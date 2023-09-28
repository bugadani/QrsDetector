//! Sliding window

use core::marker::PhantomData;

pub struct SlidingWindow<T, C>
where
    C: AsRef<[T]> + AsMut<[T]>,
{
    buffer: C,
    idx: usize,
    full: bool,
    _marker: PhantomData<T>,
}

impl<T: Default + Copy, const N: usize> Default for SlidingWindow<T, [T; N]> {
    fn default() -> Self {
        Self::new([T::default(); N])
    }
}

impl<T, C> SlidingWindow<T, C>
where
    T: Copy,
    C: AsRef<[T]> + AsMut<[T]>,
{
    pub fn new(buffer: C) -> Self {
        Self {
            buffer,
            idx: 0,
            full: false,
            _marker: PhantomData,
        }
    }

    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.idx = 0;
        self.full = false;
    }

    pub fn len(&self) -> usize {
        if self.full {
            self.buffer.as_ref().len()
        } else {
            self.idx
        }
    }

    pub fn last(&self) -> Option<T> {
        if self.idx == 0 && !self.full {
            None
        } else {
            let idx = if self.idx == 0 {
                self.buffer.as_ref().len()
            } else {
                self.idx
            };
            Some(self.buffer.as_ref()[idx - 1])
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
        let buffer = self.buffer.as_mut();
        let old = self.full.then_some(buffer[self.idx]);

        buffer[self.idx] = sample;
        self.idx = (self.idx + 1) % buffer.len();
        if self.idx == 0 {
            self.full = true;
        }

        old
    }

    pub fn iter_unordered(&self) -> impl Iterator<Item = T> + Clone + '_ {
        (0..self.len()).map(|i| self.buffer.as_ref()[i])
    }
}
