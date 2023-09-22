//! Sliding window

pub struct SlidingWindow<T, const N: usize> {
    buffer: [T; N],
    idx: usize,
    full: bool,
}

impl<T: Default + Copy, const N: usize> Default for SlidingWindow<T, N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy, const N: usize> SlidingWindow<T, N> {
    pub fn new() -> Self
    where
        T: Default,
    {
        Self {
            buffer: [T::default(); N],
            idx: 0,
            full: false,
        }
    }

    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.idx = 0;
        self.full = false;
    }

    pub fn len(&self) -> usize {
        if self.full {
            N
        } else {
            self.idx
        }
    }

    pub fn last(&self) -> Option<T> {
        if self.idx == 0 && !self.full {
            None
        } else {
            let idx = if self.idx == 0 { N } else { self.idx };
            Some(self.buffer[idx - 1])
        }
    }

    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn is_full(&self) -> bool {
        self.full
    }

    pub fn push(&mut self, sample: T) -> Option<T> {
        let old = self.full.then_some(self.buffer[self.idx]);

        self.buffer[self.idx] = sample;
        self.idx = (self.idx + 1) % self.buffer.len();
        if self.idx == 0 {
            self.full = true;
        }

        old
    }

    pub fn iter_unordered(&self) -> impl Iterator<Item = T> + Clone + '_ {
        (0..self.len()).map(|i| self.buffer[i])
    }
}
