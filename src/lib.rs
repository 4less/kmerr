#![feature(test)]
#![feature(generic_const_exprs)]

extern crate test;

pub mod kmer_utils;
pub mod utils;
pub mod consecutive;
pub mod syncmer;
pub mod minimizer;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}