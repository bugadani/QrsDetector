//! This crate provides a realtime ECG QRS detector.
//!
//! The implementation is based on [this article](https://biomedical-engineering-online.biomedcentral.com/articles/10.1186/1475-925X-3-28).
#![cfg_attr(not(test), no_std)]

#[cfg(feature = "alloc")]
extern crate alloc;

mod algorithms;
pub mod sampling;
mod sliding;

use algorithms::{F, M, R};
use sampling::SamplingFrequency;

use crate::sliding::SlidingWindow;

/// Finds QRS complexes in real-time sampled ECG signal.
///
/// # Type parameters:
///
/// - `FMW` - a buffer type containing 300ms worth of samples
/// - `FB` - a buffer type containing 50ms worth of samples
///
/// These type parameters are checked at runtime and, if incorrect, the error message will contain
/// the correct sizes.
pub struct QrsDetector<FMW, FB> {
    total_samples: u32,
    m: M,
    f: F<FMW, FB>,
    r: R,
}

impl QrsDetector<(), ()> {
    /// Creates a new QRS detector for signals sampled with `fs`. The internal buffers will be
    /// allocated on the stack as part of the `QrsDetector` structure.
    ///
    /// # Arguments
    /// * `fs` - The sampling frequency of the processed signal. For more information see
    /// [`sampling::SamplingFrequencyExt`].
    ///
    /// # Example
    /// ```rust
    /// use qrs_detector::sampling::*;
    /// use qrs_detector::QrsDetector;
    ///
    /// // Assuming 500 samples per second
    /// // Type parameters must be 300ms and 50ms in number of samples
    /// let detector = QrsDetector::new::<150, 25>(500.sps());
    /// ```
    pub fn new<const SAMPLES_300: usize, const SAMPLES_50: usize>(
        fs: SamplingFrequency,
    ) -> QrsDetector<[f32; SAMPLES_300], [f32; SAMPLES_50]> {
        QrsDetector {
            total_samples: 0,
            m: M::new(fs),
            f: F::new(fs, SlidingWindow::default(), SlidingWindow::default()),
            r: R::new(),
        }
    }

    /// Creates a new QRS detector for signals sampled with `fs`, using the provided buffers.
    ///
    /// # Arguments
    /// * `fs` - The sampling frequency of the processed signal. For more information see
    /// [`sampling::SamplingFrequencyExt`].
    /// * `f_buffer_300` - A buffer containing 300ms worth of samples.
    /// * `f_buffer_50` - A buffer containing 50ms worth of samples.
    ///
    /// # Examples
    ///
    /// The backing buffers may be slices:
    ///
    /// ```rust
    /// use qrs_detector::sampling::*;
    /// use qrs_detector::QrsDetector;
    ///
    /// // Assuming 500 samples per second
    /// // Buffers must hold 300ms and 50ms worth of samples.
    /// let mut f_buffer_300 = [0.0; 150];
    /// let mut f_buffer_50 = [0.0; 25];
    /// let detector = QrsDetector::new_from(500.sps(), &mut f_buffer_300, &mut f_buffer_50);
    /// ```
    ///
    /// The backing buffers may be arrays:
    ///
    /// ```rust
    /// use qrs_detector::sampling::*;
    /// use qrs_detector::QrsDetector;
    ///
    /// // Assuming 500 samples per second
    /// // Buffers must hold 300ms and 50ms worth of samples.
    /// let f_buffer_300 = [0.0; 150];
    /// let f_buffer_50 = [0.0; 25];
    /// let detector = QrsDetector::new_from(500.sps(), f_buffer_300, f_buffer_50);
    /// ```
    pub fn new_from<FMW, FB>(
        fs: SamplingFrequency,
        f_buffer_300: FMW,
        f_buffer_50: FB,
    ) -> QrsDetector<FMW, FB>
    where
        FMW: AsRef<[f32]> + AsMut<[f32]>,
        FB: AsRef<[f32]> + AsMut<[f32]>,
    {
        QrsDetector {
            total_samples: 0,
            m: M::new(fs),
            f: F::new(
                fs,
                SlidingWindow::new(f_buffer_300),
                SlidingWindow::new(f_buffer_50),
            ),
            r: R::new(),
        }
    }

    /// Creates a new QRS detector for signals sampled with `fs`, using the provided buffers.
    ///
    /// # Arguments
    /// * `fs` - The sampling frequency of the processed signal. For more information see
    /// [`sampling::SamplingFrequencyExt`].
    /// * `f_buffer_300` - A buffer containing 300ms worth of samples.
    /// * `f_buffer_50` - A buffer containing 50ms worth of samples.
    ///
    /// # Examples
    ///
    /// The backing buffers may be slices:
    ///
    /// ```rust
    /// use qrs_detector::sampling::*;
    /// use qrs_detector::QrsDetector;
    ///
    /// // Assuming 500 samples per second
    /// // Buffers must hold 300ms and 50ms worth of samples.
    /// let mut f_buffer_300 = [0.0; 150];
    /// let mut f_buffer_50 = [0.0; 25];
    /// let detector = QrsDetector::new_from(500.sps(), &mut f_buffer_300, &mut f_buffer_50);
    /// ```
    ///
    /// The backing buffers may be arrays:
    ///
    /// ```rust
    /// use qrs_detector::sampling::*;
    /// use qrs_detector::QrsDetector;
    ///
    /// // Assuming 500 samples per second
    /// // Buffers must hold 300ms and 50ms worth of samples.
    /// let f_buffer_300 = [0.0; 150];
    /// let f_buffer_50 = [0.0; 25];
    /// let detector = QrsDetector::new_from(500.sps(), f_buffer_300, f_buffer_50);
    /// ```
    #[cfg(feature = "alloc")]
    pub fn new_alloc(
        fs: SamplingFrequency,
    ) -> QrsDetector<alloc::boxed::Box<[f32]>, alloc::boxed::Box<[f32]>> {
        use alloc::vec;
        QrsDetector {
            total_samples: 0,
            m: M::new(fs),
            f: F::new(
                fs,
                SlidingWindow::new(vec![0.0; fs.ms_to_samples(300.0)].into_boxed_slice()),
                SlidingWindow::new(vec![0.0; fs.ms_to_samples(50.0)].into_boxed_slice()),
            ),
            r: R::new(),
        }
    }
}

impl<FMW, FB> QrsDetector<FMW, FB>
where
    FMW: AsRef<[f32]> + AsMut<[f32]>,
    FB: AsRef<[f32]> + AsMut<[f32]>,
{
    /// Resets the internal state of the detector.
    pub fn clear(&mut self) {
        self.m.clear();
        self.f.clear();
        self.r.clear();
    }

    /// Processes a sample. Returns Some sample index if a QRS complex is detected.
    pub fn update(&mut self, sample: f32) -> Option<u32> {
        self.m.update(sample);
        self.f.update(sample);
        self.r.update(self.m.current_decrement);

        let thresholds = self.thresholds();

        let result = match thresholds.total() {
            Some(mfr) if sample > mfr => {
                self.m.detection_event(sample);
                self.r.detection_event(self.total_samples);
                Some(self.total_samples)
            }
            _ => None,
        };

        self.total_samples += 1;
        result
    }

    /// Returns the current threshold value.
    /// This value is used to determine if a sample is a QRS complex.
    /// The final threshold is calculated as `M + F + R`.
    pub fn thresholds(&self) -> Thresholds {
        Thresholds {
            m: self.m.threshold(),
            f: self.f.threshold(),
            r: self.r.threshold(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Thresholds {
    pub m: Option<f32>,
    pub f: Option<f32>,
    pub r: f32,
}

impl Thresholds {
    pub fn total(&self) -> Option<f32> {
        if let (Some(m), Some(f)) = (self.m, self.f) {
            Some(m + f + self.r)
        } else {
            None
        }
    }
}
