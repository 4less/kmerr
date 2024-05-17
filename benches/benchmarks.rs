
#![feature(test)]
extern crate test;
extern crate kmerrs;

use kmerrs::*;

#[cfg(test)]
mod kmer_tests {
    use test::Bencher;

    use crate::{consecutive::kmer::{Kmer, KmerIter}, minimizer::context_free::Minimizer, syncmer::closed_syncmer::ClosedSyncmer};


    #[bench]
    fn bench_construct_kmers_iterator_option_closed_syncmer(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        const S: usize = 7; 
        let mut cs = ClosedSyncmer::<{K as u64},{S as u64}>::new();

        let mut sum = 0;
        b.iter(|| {
            let mut k_iter = KmerIter::<K>::new(&kmer.as_bytes());
            k_iter.set_assume_perfect_data(false);
            k_iter.for_each(|(_, kmer)| {
                sum += kmer.0;
                cs.load(kmer.0);
                sum += cs.is_minimizer() as u64;
            });
        })
    }
}