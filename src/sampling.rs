//! Helpers for working with sampling frequencies and sample numbers.

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct SamplingFrequency(f32);

/// Extension functions for numeric types used to create [`SamplingFrequency`] values.
///
/// # Usage
/// ```rust
/// use qrs_detector::sampling::*;
///
/// // Both values represent 500 samples per second
/// let fs = 500.sps();
/// let fs2 = 0.5.ksps();
///
/// assert_eq!(fs, fs2);
/// ```
pub trait SamplingFrequencyExt {
    fn sps(self) -> SamplingFrequency;
    fn ksps(self) -> SamplingFrequency;
}

impl SamplingFrequencyExt for f32 {
    fn sps(self) -> SamplingFrequency {
        SamplingFrequency(self)
    }

    fn ksps(self) -> SamplingFrequency {
        (self * 1000.0).sps()
    }
}

impl SamplingFrequencyExt for usize {
    fn sps(self) -> SamplingFrequency {
        SamplingFrequency(self as f32)
    }

    fn ksps(self) -> SamplingFrequency {
        (self * 1000).sps()
    }
}

impl SamplingFrequency {
    /// Convert `ms` milliseconds to number of samples
    /// ```rust
    /// # use qrs_detector::sampling::*;
    /// #
    /// let samples = 500.sps().s_to_samples(2.5);
    /// assert_eq!(samples, 1250);
    /// ```
    pub fn ms_to_samples(self, ms: f32) -> usize {
        ((ms * self.0) as usize) / 1000
    }

    /// Convert `s` seconds to number of samples
    /// ```rust
    /// # use qrs_detector::sampling::*;
    /// #
    /// let samples = 500.sps().ms_to_samples(500.0);
    /// assert_eq!(samples, 250);
    /// ```
    pub fn s_to_samples(self, s: f32) -> usize {
        self.ms_to_samples(s * 1000.0)
    }
}
