pub trait Minimizer {
    fn is_minimizer(&mut self, hash: u64) -> bool;
}