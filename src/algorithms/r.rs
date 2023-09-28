use crate::sliding::SlidingWindow;

#[derive(Copy, Clone, Debug)]
enum RState {
    Ignore,
    InitBuffer,
    NoDecrease(u32, u32),    // samples remaining, average rr interval
    Decrease(u32, f32, f32), // samples remaining, value, decrement
    Constant(f32),
}

pub struct R {
    state: RState,
    rr: SlidingWindow<u32, [u32; 5]>,
    prev_idx: u32, // no need to make it an Option
}

impl R {
    pub fn new() -> Self {
        Self {
            state: RState::Ignore,
            rr: SlidingWindow::default(),
            prev_idx: 0,
        }
    }

    fn enter_no_decrease(&mut self) {
        let rr_sum: u32 = self.rr.iter_unordered().sum();
        let rr_avg = rr_sum / 5;
        self.state = RState::NoDecrease(rr_avg * 2 / 3, rr_avg);
    }

    pub fn update(&mut self, m_decrement: f32) {
        self.state = match self.state {
            RState::NoDecrease(0, rr_avg) => RState::Decrease(rr_avg / 3, 0.0, m_decrement / 1.4),
            RState::NoDecrease(samples, rr_avg) => RState::NoDecrease(samples - 1, rr_avg),
            RState::Decrease(0, r, _) => RState::Constant(r),
            RState::Decrease(samples, r, decrement) => {
                RState::Decrease(samples - 1, r - decrement, decrement)
            }
            o => o,
        };
    }

    pub fn threshold(&self) -> f32 {
        match self.state {
            RState::Ignore | RState::InitBuffer | RState::NoDecrease(_, _) => 0.0,
            RState::Constant(r) | RState::Decrease(_, r, _) => r,
        }
    }

    pub fn detection_event(&mut self, idx: u32) {
        match self.state {
            RState::Ignore => self.state = RState::InitBuffer,
            RState::InitBuffer => {
                self.rr.push(idx.wrapping_sub(self.prev_idx));
                if self.rr.is_full() {
                    self.enter_no_decrease();
                }
            }
            _ => {
                self.rr.push(idx.wrapping_sub(self.prev_idx));
                self.enter_no_decrease();
            }
        };

        self.prev_idx = idx;
    }
}
