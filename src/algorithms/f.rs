use crate::{sampling::*, sliding::SlidingWindow};

#[derive(Copy, Clone, Debug)]
enum FState {
    Ignore(usize),
    Init(usize, f32),
    Integrate(f32),
}

pub struct F<const SAMPLES_350: usize, const SAMPLES_50: usize> {
    /// F should be initialized at the same time as M is, skip earlier samples
    state: FState,
    /// 350ms of the individual max samples of the 50ms buffer
    f_max_window: SlidingWindow<f32, SAMPLES_350>,
    /// 50ms window of the signal
    f_buffer: SlidingWindow<f32, SAMPLES_50>,
}

impl<const SAMPLES_350: usize, const SAMPLES_50: usize> F<SAMPLES_350, SAMPLES_50> {
    pub fn new(fs: SamplingFrequency) -> Self {
        // sanity check buffer sizes
        debug_assert_eq!(
            SAMPLES_350,
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
            f_max_window: SlidingWindow::new(),
            f_buffer: SlidingWindow::new(),
        }
    }

    fn update_f_buffers(&mut self, sample: f32) -> (Option<f32>, f32) {
        // TODO: there are some special cases where the max search can be skipped
        self.f_buffer.push(sample);

        // Calculate maximum value in the latest 50ms window
        let max = self
            .f_buffer
            .iter_unordered()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        // Keep the 50ms maximum values for each sample in latest 300ms window
        // The oldest sample corresponds to the oldest 50ms in the latest 350ms window
        // TODO FIXME: off by some error :)
        let old = self.f_max_window.push(max);

        (old, max)
    }

    pub fn update(&mut self, sample: f32) {
        self.state = match self.state {
            FState::Ignore(1) => FState::Init(SAMPLES_350 - 1, 0.0),
            FState::Ignore(n) => FState::Ignore(n - 1),
            FState::Init(n, favg) => {
                let favg = favg + sample;
                self.update_f_buffers(sample);

                if n == 0 {
                    FState::Integrate(favg / (SAMPLES_350 as f32))
                } else {
                    FState::Init(n - 1, favg)
                }
            }
            FState::Integrate(f) => {
                let (oldest_max, max) = self.update_f_buffers(sample);
                FState::Integrate(f + (max - oldest_max.unwrap()) / 150.0)
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
