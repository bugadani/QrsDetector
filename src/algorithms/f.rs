use crate::{sampling::*, sliding::SlidingWindow};

#[derive(Copy, Clone, Debug)]
enum FState {
    Ignore(usize),
    Init(usize, f32),
    Integrate(f32),
}

pub struct F<FMW, FB> {
    fs: SamplingFrequency,

    /// F should be initialized at the same time as M is, skip earlier samples
    state: FState,

    /// 300ms of the individual max samples of the 50ms buffer.
    /// We store the maxima to save some memory and computation.
    f_max_window: SlidingWindow<f32, FMW>,

    /// 50ms window of the signal
    f_buffer: SlidingWindow<f32, FB>,
}

impl<FMW, FB> F<FMW, FB>
where
    FMW: AsRef<[f32]> + AsMut<[f32]>,
    FB: AsRef<[f32]> + AsMut<[f32]>,
{
    pub fn new(
        fs: SamplingFrequency,
        f_max_window: SlidingWindow<f32, FMW>,
        f_buffer: SlidingWindow<f32, FB>,
    ) -> Self {
        // sanity check buffer sizes
        debug_assert_eq!(
            f_max_window.capacity(),
            fs.ms_to_samples(300.0),
            "Incorrect FMW type parameter. Buffer must be {} samples long.",
            fs.ms_to_samples(300.0)
        );
        debug_assert_eq!(
            f_buffer.capacity(),
            fs.ms_to_samples(50.0),
            "Incorrect FB type parameter. Buffer must be {} samples long.",
            fs.ms_to_samples(50.0)
        );

        Self {
            fs,
            state: FState::Ignore(fs.s_to_samples(2.65)),
            f_max_window,
            f_buffer,
        }
    }

    pub fn clear(&mut self) {
        self.state = FState::Ignore(self.fs.s_to_samples(2.65));
        self.f_max_window.clear();
        self.f_buffer.clear();
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
            FState::Ignore(1) => FState::Init(self.f_max_window.capacity() - 1, 0.0),
            FState::Ignore(n) => FState::Ignore(n - 1),
            FState::Init(n, favg) => {
                let favg = favg + sample;
                self.update_f_buffers(sample);

                if n == 0 {
                    FState::Integrate(favg.max(0.0) / (self.f_max_window.capacity() as f32))
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
