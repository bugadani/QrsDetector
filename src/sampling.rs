/// Type alias for clarity. For more information see [`SamplingFreqencyExt`](./trait.SamplingFrequencyExt.html)
/// and [`SamplingFrequencyHelpers`](./trait.SamplingFrequencyHelpers.html).
pub type SamplingFrequency = f32;

/// Extension functions for numeric types used to create [`SamplingFreqency`](./struct.SamplingFrequency.html) values.
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
        self
    }

    fn ksps(self) -> SamplingFrequency {
        (self * 1000.0).sps()
    }
}

impl SamplingFrequencyExt for usize {
    fn sps(self) -> SamplingFrequency {
        self as SamplingFrequency
    }

    fn ksps(self) -> SamplingFrequency {
        (self * 1000).sps()
    }
}

/// Helper functions to make some sampling time related conversions simpler.
pub trait SamplingFrequencyHelpers {
    /// Convert `ms` milliseconds to number of samples
    /// ```rust
    /// # use qrs_detector::sampling::*;
    /// #
    /// let samples = 500.sps().s_to_samples(2.5);
    /// assert_eq!(samples, 1250);
    /// ```
    fn ms_to_samples(self, ms: f32) -> u32;

    /// Convert `s` seconds to number of samples
    /// ```rust
    /// # use qrs_detector::sampling::*;
    /// #
    /// let samples = 500.sps().ms_to_samples(500.0);
    /// assert_eq!(samples, 250);
    /// ```
    fn s_to_samples(self, s: f32) -> u32;
}

impl SamplingFrequencyHelpers for SamplingFrequency {
    fn ms_to_samples(self, ms: f32) -> u32 {
        ((ms * self) as u32) / 1000
    }
    fn s_to_samples(self, s: f32) -> u32 {
        self.ms_to_samples(s * 1000.0)
    }
}
