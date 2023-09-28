use crate::{sampling::*, sliding::SlidingWindow};

#[derive(Copy, Clone, Debug)]
enum FState {
    Ignore(usize),
    Init(usize, f32),
    Integrate(f32),
}

pub struct F<const SAMPLES_300: usize, const SAMPLES_50: usize> {
    /// F should be initialized at the same time as M is, skip earlier samples
    state: FState,
    /// 300ms of the individual max samples of the 50ms buffer.
    /// We store the maxima to save some memory and computation.
    f_max_window: SlidingWindow<f32, [f32; SAMPLES_300]>,
    /// 50ms window of the signal
    f_buffer: SlidingWindow<f32, [f32; SAMPLES_50]>,
}

impl<const SAMPLES_300: usize, const SAMPLES_50: usize> F<SAMPLES_300, SAMPLES_50> {
    pub fn new(fs: SamplingFrequency) -> Self {
        // sanity check buffer sizes
        debug_assert_eq!(
            SAMPLES_300,
            fs.ms_to_samples(300.0),
            "Incorrect type parameters, must be <{}, {}>",
            fs.ms_to_samples(300.0),
            fs.ms_to_samples(50.0)
        );
        debug_assert_eq!(
            SAMPLES_50,
            fs.ms_to_samples(50.0),
            "Incorrect type parameters, must be <{}, {}>",
            fs.ms_to_samples(300.0),
            fs.ms_to_samples(50.0)
        );

        Self {
            state: FState::Ignore(fs.s_to_samples(2.65)),
            f_max_window: SlidingWindow::default(),
            f_buffer: SlidingWindow::default(),
        }
    }

    fn update_f_buffers(&mut self, sample: f32) -> (Option<f32>, f32) {
        // TODO: there are some special cases where the max search can be skipped
        self.f_buffer.push(sample);

        // Calculate maximum value in the latest 50ms window
        let max: f32 = self
            .f_buffer
            .iter_unordered()
            .fold(0.0, |acc, x| acc.max(x));

        // Keep the 50ms maximum values for each sample in latest 300ms window
        // The oldest sample corresponds to the oldest 50ms in the latest 350ms window
        // TODO FIXME: off by some error :)
        let old = self.f_max_window.push(max);

        (old, max)
    }

    pub fn update(&mut self, sample: f32) {
        self.state = match self.state {
            FState::Ignore(1) => FState::Init(SAMPLES_300 - 1, 0.0),
            FState::Ignore(n) => FState::Ignore(n - 1),
            FState::Init(n, favg) => {
                let favg = favg + sample;
                self.update_f_buffers(sample);

                if n == 0 {
                    FState::Integrate(favg.max(0.0) / (SAMPLES_300 as f32))
                } else {
                    FState::Init(n - 1, favg)
                }
            }
            FState::Integrate(f) => {
                let (oldest_max, max) = self.update_f_buffers(sample);
                FState::Integrate((f + (max - oldest_max.unwrap()) / 150.0).max(0.0))
            }
        };
    }

    pub fn threshold(&self) -> Option<f32> {
        if let FState::Integrate(f) = self.state {
            Some(f)
        } else {
            None
        }
    }
}
