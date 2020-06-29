//! This crate provides a realtime ECG QRS detector.
//!
//! The implementation is based on [this article](https://biomedical-engineering-online.biomedcentral.com/articles/10.1186/1475-925X-3-28).
//!
//! `QrsDetector` does not use dynamic memory allocation. Instead, this crate relies on
//! [`generic_array`](../generic_array/index.html),
//! which means there's a bit of setup necessary to set the size of internal buffers.
//!
//! For more information, see [`QrsDetector`](./struct.QrsDetector.html).
#![cfg_attr(not(test), no_std)]

mod internal {
    use micromath::F32Ext;

    use sliding_window::{typenum::consts::U2, SlidingWindow};

    pub fn max(a: f32, b: f32) -> f32 {
        if a > b {
            a
        } else {
            b
        }
    }

    pub struct Differentiator {
        window: SlidingWindow<f32, U2>,
    }

    impl Differentiator {
        pub fn new() -> Self {
            Self {
                window: SlidingWindow::new(),
            }
        }

        pub fn update(&mut self, sample: f32) -> Option<f32> {
            self.window.insert(sample).map(|old| (sample - old).abs())
        }
    }

    #[cfg(test)]
    mod test_diff {
        use super::Differentiator;

        #[test]
        fn sanity_check() {
            let mut diff = Differentiator::new();

            assert_eq!(None, diff.update(1.0));
            assert_eq!(None, diff.update(2.0));
            assert_eq!(Some(2.0), diff.update(3.0));
            assert_eq!(Some(2.0), diff.update(4.0));
            assert_eq!(Some(0.0), diff.update(3.0));
            assert_eq!(Some(1.0), diff.update(3.0));
            assert_eq!(Some(1.0), diff.update(2.0));
            assert_eq!(Some(2.0), diff.update(1.0));
        }
    }
}

use if_chain::if_chain;
use sliding_window::Size;

use internal::Differentiator;
pub use sliding_window::typenum;

mod algorithms {
    use sliding_window::{typenum::consts::U5, Size, SlidingWindow};

    use crate::internal::max;
    use crate::sampling::*;

    #[derive(Copy, Clone, Debug)]
    pub enum MState {
        Init(u32, f32),
        Disallow(u32, f32),
        Decreasing(u32, f32, f32),
        ConstantLow(f32),
    }

    pub struct M {
        state: MState,
        mm: SlidingWindow<f32, U5>,
        fs: SamplingFrequency,
        pub current_decrement: f32,
    }

    impl M {
        pub fn new(fs: SamplingFrequency) -> Self {
            Self {
                fs,
                mm: SlidingWindow::new(),
                state: MState::Init(fs.s_to_samples(3.0), 0.0),
                current_decrement: 0.0,
            }
        }

        fn m(&self) -> f32 {
            self.mm.iter_unordered().sum()
        }

        pub fn update(&mut self, sample: f32) -> Option<f32> {
            self.state = match self.state {
                MState::Init(0, m) => {
                    let m = 0.6 * max(m, sample) / 5.0; // divide by 5 for averaging

                    for _ in 0..5 {
                        self.mm.insert(m);
                    }

                    let n_samples = self.fs.s_to_samples(1.0);
                    let decrement = m * 0.4 / n_samples as f32;
                    self.current_decrement = decrement;
                    MState::Decreasing(n_samples, m, decrement)
                }
                MState::Init(samples, m) => MState::Init(samples - 1, max(m, sample)),
                MState::Disallow(0, m) => {
                    let m = 0.6 * max(m, sample) / 5.0; // divide by 5 for averaging
                    let prev_m = self.mm[4];
                    let new_m = if m > prev_m * 1.5 { 1.1 * prev_m } else { m };
                    self.mm.insert(new_m);

                    let m = self.m();
                    let n_samples = self.fs.s_to_samples(1.0);
                    let decrement = m * 0.4 / n_samples as f32;
                    self.current_decrement = decrement;
                    MState::Decreasing(n_samples, m, decrement)
                }
                MState::Disallow(samples, m) => MState::Disallow(samples - 1, max(m, sample)),
                MState::Decreasing(0, m, _) => MState::ConstantLow(m),
                MState::Decreasing(samples, m, decrease_amount) => {
                    MState::Decreasing(samples - 1, m - decrease_amount, decrease_amount)
                }
                MState::ConstantLow(m) => MState::ConstantLow(m),
            };

            match self.state {
                MState::Init(_, _) | MState::Disallow(_, _) => None,
                MState::Decreasing(_, m, _) | MState::ConstantLow(m) => Some(m),
            }
        }

        pub fn detection_event(&mut self, sample: f32) {
            self.state = MState::Disallow(self.fs.s_to_samples(0.225), sample);
        }
    }

    #[derive(Copy, Clone, Debug)]
    pub enum FState {
        Ignore(u32),
        Init(u32, f32),
        Integrate(f32),
    }

    pub struct F<FMW, FB>
    where
        FMW: Size<f32>,
        FB: Size<f32>,
    {
        /// Sampling frequency
        fs: SamplingFrequency,
        /// F should be initialized at the same time as M is, skip earlier samples
        state: FState,
        /// 350ms of the individual max samples of the 50ms buffer
        f_max_window: SlidingWindow<f32, FMW>,
        /// 50ms window of the signal
        f_buffer: SlidingWindow<f32, FB>,
    }

    impl<FMW, FB> F<FMW, FB>
    where
        FMW: Size<f32>,
        FB: Size<f32>,
    {
        pub fn new(fs: SamplingFrequency) -> Self {
            // sanity check buffer sizes
            debug_assert_eq!(
                FMW::to_u32(),
                fs.ms_to_samples(350.0),
                "Incorrect type parameters, must be <U{}, U{}>",
                fs.ms_to_samples(350.0),
                fs.ms_to_samples(50.0)
            );
            debug_assert_eq!(
                FB::to_u32(),
                fs.ms_to_samples(50.0),
                "Incorrect type parameters, must be <U{}, U{}>",
                fs.ms_to_samples(350.0),
                fs.ms_to_samples(50.0)
            );

            Self {
                fs,
                state: FState::Ignore(fs.s_to_samples(2.65)),
                f_max_window: SlidingWindow::new(),
                f_buffer: SlidingWindow::new(),
            }
        }

        fn update_f_buffers(&mut self, sample: f32) -> (Option<f32>, f32) {
            // TODO: there are some special cases where the max search can be skipped
            self.f_buffer.insert(sample);
            let max = *self
                .f_buffer
                .iter_unordered()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            let old = self.f_max_window.insert(max);

            (old, max)
        }

        pub fn update(&mut self, sample: f32) -> Option<f32> {
            self.state = match self.state {
                FState::Ignore(1) => FState::Init(self.fs.ms_to_samples(350.0), 0.0),
                FState::Ignore(n) => FState::Ignore(n - 1),
                FState::Init(1, favg) => {
                    let favg = favg + sample;
                    self.update_f_buffers(sample);

                    FState::Integrate(favg / (FMW::to_u32() as f32))
                }
                FState::Init(n, favg) => {
                    let favg = favg + sample;
                    self.update_f_buffers(sample);

                    FState::Init(n - 1, favg)
                }
                FState::Integrate(f) => {
                    let (oldest_max, max) = self.update_f_buffers(sample);
                    FState::Integrate(f + (max - oldest_max.unwrap()) / 150.0)
                }
            };

            if let FState::Integrate(f) = self.state {
                Some(f)
            } else {
                None
            }
        }
    }

    #[derive(Copy, Clone, Debug)]
    pub enum RState {
        Ignore,
        InitBuffer,
        NoDecrease(u32, u32),    // samples remaining, average rr interval
        Decrease(u32, f32, f32), // samples remaining, value, decrement
        Constant(f32),
    }

    pub struct R {
        state: RState,
        rr: SlidingWindow<u32, U5>,
        prev_idx: u32, // no need to make it an Option
    }

    impl R {
        pub fn new() -> Self {
            Self {
                state: RState::Ignore,
                rr: SlidingWindow::new(),
                prev_idx: 0,
            }
        }

        fn enter_no_decrease(&mut self) {
            let rr_sum: u32 = self.rr.iter().sum();
            let rr_avg = rr_sum / 5;
            self.state = RState::NoDecrease(rr_avg * 2 / 3, rr_avg);
        }

        pub fn update(&mut self, m_decrement: f32) -> f32 {
            self.state = match self.state {
                RState::NoDecrease(0, rr_avg) => {
                    RState::Decrease(rr_avg / 3, 0.0, m_decrement / 1.4)
                }
                RState::NoDecrease(samples, rr_avg) => RState::NoDecrease(samples - 1, rr_avg),
                RState::Decrease(0, r, _) => RState::Constant(r),
                RState::Decrease(samples, r, decrement) => {
                    RState::Decrease(samples - 1, r - decrement, decrement)
                }
                o => o,
            };

            match self.state {
                RState::Ignore | RState::InitBuffer | RState::NoDecrease(_, _) => 0.0,
                RState::Constant(r) | RState::Decrease(_, r, _) => r,
            }
        }

        pub fn detection_event(&mut self, idx: u32) {
            match self.state {
                RState::Ignore => self.state = RState::InitBuffer,
                RState::InitBuffer => {
                    self.rr.insert(idx.wrapping_sub(self.prev_idx));
                    if self.rr.is_full() {
                        self.enter_no_decrease();
                    }
                }
                _ => {
                    self.rr.insert(idx.wrapping_sub(self.prev_idx));
                    self.enter_no_decrease();
                }
            };

            self.prev_idx = idx;
        }
    }
}

/// Helpers for working with sampling frequencies and sample numbers.
pub mod sampling;
use sampling::SamplingFrequency;

use algorithms::{F, M, R};

/// Find QRS complex in real-time sampled ECG signal.
///
/// # Type parameters:
///
/// The `QrsDetector` is built upon [`generic_array`](../generic_array/index.html).
/// This has an unfortunate implementation detail where the internal buffer sizes must be set on
/// the type level.
///
/// - `FMW` - number of samples representing 350ms, as [`typenum::U*`](../typenum/consts/index.html)
/// - `FB` - number of samples representing 50ms, as [`typenum::U*`](../typenum/consts/index.html)
///
/// These type parameters are checked at runtime and if incorrect and the error message will contain
/// the correct sizes.
pub struct QrsDetector<FMW, FB>
where
    FMW: Size<f32>,
    FB: Size<f32>,
{
    fs: SamplingFrequency,
    differentiator: Differentiator,
    total_samples: u32,
    m: M,
    f: F<FMW, FB>,
    r: R,
}

impl<FMW, FB> QrsDetector<FMW, FB>
where
    FMW: Size<f32>,
    FB: Size<f32>,
{
    /// Returns a new QRS detector for signals sampled at `fs` sampling frequency.
    ///
    /// # Arguments
    /// * `fs` - The sampling frequency of the processed signal. For more information see
    /// [`SamplingFrequencyExt`](./sampling/trait.SamplingFrequencyExt.html).
    ///
    /// # Example
    /// ```rust
    /// use qrs_detector::sampling::*;
    /// use qrs_detector::QrsDetector;
    /// use qrs_detector::typenum::{U175, U25};
    ///
    /// // Assuming 500 samples per second
    /// // Type parameters must be 350ms and 50ms in number of samples
    /// let detector: QrsDetector<U175, U25> = QrsDetector::new(500.sps());
    /// ```
    pub fn new(fs: SamplingFrequency) -> Self {
        Self {
            fs,
            total_samples: 0,
            differentiator: Differentiator::new(),
            m: M::new(fs),
            f: F::new(fs),
            r: R::new(),
        }
    }

    /// Reset the internal state of the detector
    pub fn clear(&mut self) {
        *self = Self::new(self.fs);
    }

    /// Process a sample. Returns Some sample index if a QRS complex is detected.
    pub fn update(&mut self, input: f32) -> Option<u32> {
        let result = self.differentiator.update(input).and_then(|sample| {
            let m = self.m.update(sample);
            let f = self.f.update(sample);
            let r = self.r.update(self.m.current_decrement);

            if_chain! {
                if let Some(m) = m;
                if let Some(f) = f;
                if sample > m + f + r;
                then {
                    self.m.detection_event(sample);
                    self.r.detection_event(self.total_samples);
                    Some(self.total_samples)
                } else {
                    None
                }
            }
        });

        self.total_samples += 1;

        result
    }
}
