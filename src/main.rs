use std::{cmp::min, hash::{DefaultHasher, Hash, Hasher}};

use kmerrs::{consecutive::kmer::KmerIter, minimizer::{context_free::Minimizer, universal_minimizer::UniMin}, syncmer::closed_syncmer::ClosedSyncmer};
use kmers::KmerIterator;

fn main() {
    let kmer: &str = "GACGCTCGACTGGTCTGCCZGTAGCCGTACGACTCTCTGTCAATGTTTGCTTCTCTTCATTATTATATTGATCTACGCTACTACTAGCTGATTATATATATAGACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        

    const K: usize = 15;

    let kit = KmerIterator::new(K, kmer.as_bytes().iter());
    println!("{}", kmer);
    for (pos, (kmer1, kmer2)) in kit.enumerate() {
        println!("{}{} rev{}", " ".repeat(pos), kmer1.render(K), kmer2.render(K));
    }

    let mut state = DefaultHasher::new();
    println!("\n\n\n");
    // let mut cs = ClosedSyncmer::<15, 4>::new();
    let mut cs = UniMin::<15, 8>::new();
    let mut kit = KmerIter::<K, true>::new(kmer.as_bytes());


    println!("{}", kmer);
    let mut smer_count = 0;
    for (pos, fwd, rev) in kit {
        let kmer = min(fwd, rev);
        kmer.hash(&mut state);
        let hash = state.finish();
        if !cs.is_minimizer(hash >> (64-(K*2))) {
            continue
        }
        smer_count += 1;
        println!("{}{}", " ".repeat(pos), kmer.to_string().unwrap());
    }
    println!("Compression rate: {}", smer_count as f64/(kmer.len() - K + 1) as f64);
}