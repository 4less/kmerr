use super::context_free::Minimizer;

pub struct UniMin<const K: usize, const C: usize> {
    threshold: u64,
}

impl<const K: usize, const C: usize> UniMin<K, C> {
    pub fn new() -> Self {
        Self {
            threshold: dbg!(((1 << (K * 2)) -1)/C as u64),
        }
    }
}

impl<const K: usize, const C: usize> Minimizer for UniMin<K, C> {
    fn is_minimizer(&mut self, hash: u64) -> bool {
        hash < self.threshold
    }
}