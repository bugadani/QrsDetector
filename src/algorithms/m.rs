use crate::{sampling::*, sliding::SlidingWindow};

#[derive(Copy, Clone, Debug)]
enum MState {
    Init(usize, f32),
    Disallow(usize, f32),
    Decreasing(usize, f32, f32),
    ConstantLow(f32),
}

pub struct M {
    state: MState,
    mm: SlidingWindow<f32, 5>,
    fs: SamplingFrequency,
    pub current_decrement: f32,
}

impl M {
    pub fn new(fs: SamplingFrequency) -> Self {
        Self {
            fs,
            mm: SlidingWindow::new(),
            // Initially M = 0.6*max(Y) is set for the first 3 s [originally 5s] of the signal
            state: MState::Init(fs.s_to_samples(3.0), 0.0),
            current_decrement: 0.0,
        }
    }

    fn m(&self) -> f32 {
        // M is calculated as an average value of MM.
        // Divide by 5 was done while calculating the individual Mx values
        self.mm.iter_unordered().sum()
    }

    pub fn update(&mut self, sample: f32) {
        self.state = match self.state {
            MState::Init(0, m) => {
                // Initially M = 0.6*max(Y) is set for the first 3 s [originally 5s] of the signal
                // A buffer with 5 steep-slope threshold values is preset:
                // MM = [M1 M2 M3 M4 M5],
                // where M1 ÷ M5 are equal to M
                let m = 0.6 * m.max(sample);

                for _ in 0..5 {
                    self.mm.push(m / 5.0);
                }

                // It is not clear in the article what to do initially:
                // - wait for a QRS detection and then start the M algorithm with Disallow state
                // - decrease using "low slope" immediately
                // This implementation uses the second option

                // M is decreased in an interval 225 to 1225 ms [originally 200 to 1200 ms]
                // following the last QRS detection at a low slope, reaching 60 % of its
                // refreshed value at 1225 ms [originally 1200 ms].
                let n_samples = self.fs.s_to_samples(1.0);
                let decrement = m * 0.4 / n_samples as f32;
                self.current_decrement = decrement;

                MState::Decreasing(n_samples, m, decrement)
            }

            // Initially M = 0.6*max(Y) is set for the first 3 s [originally 5s] of the signal
            // Collect maximum value while in Init state
            MState::Init(samples, m) => MState::Init(samples - 1, m.max(sample)),

            MState::Disallow(0, m) => {
                // In the interval QRS ÷ QRS+200ms a new value of M5 is calculated:
                // newM 5 = 0.6*max(Yi)
                let m = 0.6 * m.max(sample) / 5.0; // divide by 5 for averaging

                // The estimated newM 5 value can become quite high, if steep slope premature
                // ventricular contraction or artifact appeared, and for that reason it is
                // limited to newM5 = 1.1* M5 if newM 5 > 1.5* M5.
                let prev_m = self.mm.last().unwrap_or(0.0);
                if m > prev_m * 1.5 {
                    self.mm.push(1.1 * prev_m);
                } else {
                    self.mm.push(m);
                }

                let m = self.m();

                // M is decreased in an interval 225 to 1225 ms [originally 200 to 1200 ms]
                // following the last QRS detection at a low slope, reaching 60 % of its
                // refreshed value at 1225 ms [originally 1200 ms].
                let n_samples = self.fs.s_to_samples(1.0);
                let decrement = m * 0.4 / n_samples as f32;
                self.current_decrement = decrement;

                MState::Decreasing(n_samples, m, decrement)
            }

            // In the interval QRS ÷ QRS+200ms a new value of M5 is calculated:
            // newM 5 = 0.6*max(Yi)
            // Collect maximum value while in Disallow state
            MState::Disallow(samples, m) if sample > m => {
                // if we found a new maximum, extend the disallow period
                MState::Disallow(samples.max(self.fs.s_to_samples(0.2)), sample)
            }
            MState::Disallow(samples, m) => MState::Disallow(samples - 1, m),

            // After 1225 ms [originally 1200 ms] M remains unchanged.
            MState::Decreasing(0, m, _) => MState::ConstantLow(m),

            // M is decreased in an interval 225 to 1225 ms [originally 200 to 1200 ms]
            // following the last QRS detection at a low slope, reaching 60 % of its
            // refreshed value at 1225 ms [originally 1200 ms].
            MState::Decreasing(samples, m, decrease_amount) => {
                // Linear decrease using precomputed decrement value
                MState::Decreasing(samples - 1, m - decrease_amount, decrease_amount)
            }

            // After 1225 ms [originally 1200 ms] M remains unchanged.
            MState::ConstantLow(m) => MState::ConstantLow(m),
        };
    }

    pub fn threshold(&self) -> Option<f32> {
        match self.state {
            MState::Init(_, _) | MState::Disallow(_, _) => None,
            MState::Decreasing(_, m, _) | MState::ConstantLow(m) => Some(m),
        }
    }

    pub fn detection_event(&mut self, sample: f32) {
        // No detection is allowed 225 ms [originally 200 ms] after the current one.
        self.state = MState::Disallow(self.fs.s_to_samples(0.225), sample);
    }
}
