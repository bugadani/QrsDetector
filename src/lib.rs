//! This crate provides a realtime ECG QRS detector.
//!
//! Implementation based on https://biomedical-engineering-online.biomedcentral.com/articles/10.1186/1475-925X-3-28
//!
//! Because this crate relies on generic_array, there's a bit of setup necessary to set the size of internal buffers.
#![cfg_attr(not(test), no_std)]

mod internal {
    use micromath::F32Ext;

    use sliding_window::{
        SlidingWindow, Producer,
        typenum::consts::U2
    };

    pub struct Differentiator {
        window: SlidingWindow<f32, U2>
    }

    impl Differentiator {
        pub fn new() -> Self {
            Self {
                window: SlidingWindow::new()
            }
        }

        pub fn clear(&mut self) {
            self.window.clear();
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

            diff.clear();
            assert_eq!(None, diff.update(1.0));
        }
    }
}

mod sampling {
    pub type SamplingFrequency = f32;

    pub trait SamplingFrequencyExt {
        fn sps(self) -> SamplingFrequency;
        fn ksps(self) -> SamplingFrequency;
    }

    pub trait SamplingFrequencyHelpers {
        fn ms_to_samples(self, ms: f32) -> u32;
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

    impl SamplingFrequencyHelpers for SamplingFrequency {
        fn ms_to_samples(self, ms: f32) -> u32 {
            ((ms * self) as u32) / 1000
        }
    }

    #[cfg(test)]
    mod sps_test {
        use super::SamplingFrequencyExt;
        use super::SamplingFrequencyHelpers;

        #[test]
        fn sanity_test() {
            assert_eq!(50, 1.ksps().ms_to_samples(50.0));
            assert_eq!(50, 1.0.ksps().ms_to_samples(50.0));

            assert_eq!(1, 100.sps().ms_to_samples(10.0));
        }
    }
}