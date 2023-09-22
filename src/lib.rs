//! This crate provides a realtime ECG QRS detector.
//!
//! The implementation is based on [this article](https://biomedical-engineering-online.biomedcentral.com/articles/10.1186/1475-925X-3-28).
#![cfg_attr(not(test), no_std)]

mod algorithms;
pub mod sampling;
mod sliding;

use algorithms::{F, M, R};
use if_chain::if_chain;
use sampling::SamplingFrequency;

/// Finds QRS complexes in real-time sampled ECG signal.
///
/// # Type parameters:
///
/// The internal buffer sizes must be set on the type level.
///
/// - `FMW` - number of samples representing 350ms
/// - `FB` - number of samples representing 50ms
///
/// These type parameters are checked at runtime and if incorrect and the error message will contain
/// the correct sizes.
pub struct QrsDetector<const SAMPLES_350: usize, const SAMPLES_50: usize> {
    fs: SamplingFrequency,
    total_samples: u32,
    m: M,
    f: F<SAMPLES_350, SAMPLES_50>,
    r: R,
}

impl<const SAMPLES_350: usize, const SAMPLES_50: usize> QrsDetector<SAMPLES_350, SAMPLES_50> {
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
    ///
    /// // Assuming 500 samples per second
    /// // Type parameters must be 300ms and 50ms in number of samples
    /// let detector: QrsDetector<150, 25> = QrsDetector::new(500.sps());
    /// ```
    pub fn new(fs: SamplingFrequency) -> Self {
        Self {
            fs,
            total_samples: 0,
            m: M::new(fs),
            f: F::new(fs),
            r: R::new(),
        }
    }

    /// Resets the internal state of the detector.
    pub fn clear(&mut self) {
        *self = Self::new(self.fs);
    }

    /// Processes a sample. Returns Some sample index if a QRS complex is detected.
    pub fn update(&mut self, sample: f32) -> Option<u32> {
        let m = self.m.update(sample);
        let f = self.f.update(sample);
        let r = self.r.update(self.m.current_decrement);

        let result = if_chain! {
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
        };

        self.total_samples += 1;
        result
    }
}
