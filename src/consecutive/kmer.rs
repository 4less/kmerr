
use std::cmp::min;

use crate::kmer_utils::{allowed, canonical_rc, canonical_rc2, canonical_to_u8, next_imperfect_nuc, nuc_perfect, slice_to_canonical_4, NUC_U8_TO_BIT, NUC_U8_TO_BIT_RC};
use kmers::{Kmer as ExtKmer};

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Debug)]
pub struct Kmer<const K: usize>(pub u64);

impl<const K: usize> Kmer<K> {
    #[inline]
    pub fn from_slice(slice: &[u8]) -> Option<Kmer<K>> {
        assert!(slice.len() == K);
        Some(Kmer(slice_to_canonical_4(slice)))
    }

    #[inline]
    pub fn middle<const M: usize>(&self) -> Kmer<M> {
        assert!(M < K);
        assert_eq!(M & 0b1, K & 0b1);
        let mask: u64 = (1 << M*2) - 1;
        let shift = K - M;
        Kmer::<M>((self.0 & (mask << shift)) >> shift)
    }

    #[inline]
    pub fn flanks<const F: usize>(&self) -> Kmer<F> {
        assert!(F < K);
        assert_eq!(F & 0b1, 0);

        let mask: u64 = (1 << F) - 1;
        let right = self.0 & mask;
        let left = self.0 & (mask << (2*K - F));
        let fmer = right | (left >> (2*(K - F)));

        Kmer::<F>(fmer)
    }

    #[inline]
    pub fn to_string(&self) -> Result<String, std::string::FromUtf8Error> {
        String::from_utf8(canonical_to_u8::<K>(self.0).to_vec())
    }   

    #[inline]
    pub fn rc(&self) -> Kmer<K> {
        Kmer(canonical_rc2::<K>(self.0))
    }


    #[inline]
    pub fn is_smallest_rc(&self) -> bool {
        self.0 == std::cmp::min(self.0, self.rc().0)
    }

    #[inline]
    pub fn get_smallest_rc(&self) -> Kmer<K> {
        if self.is_smallest_rc() {
            self.clone()
        } else {
            self.rc()
        }
    }

    #[inline]
    pub fn append_option(&self, nuc: u8) -> Option<Kmer<K>> {
        // if !allowed(nuc) { return None; };
        Some(Kmer(( (self.0 << 2) | NUC_U8_TO_BIT[nuc as usize - 65] as u64) & ((1u64 << K*2) - 1)))
    }

    #[inline]
    pub fn prepend_rc_option(&self, nuc: u8) -> Option<Kmer<K>> {
        // if !allowed(nuc) { return None; };
        Some(Kmer(self.0 >> 2 | ((NUC_U8_TO_BIT_RC[nuc as usize - 65] as u64) << (K-1)*2)))
    }

    #[inline]
    pub fn append(&self, nuc: u8) -> Kmer<K> {
        Kmer(( (self.0 << 2) | NUC_U8_TO_BIT[nuc as usize - 65] as u64) & ((1u64 << K*2) - 1))
    }

    #[inline]
    pub fn prepend_rc(&self, nuc: u8) -> Kmer<K> {
        Kmer(self.0 >> 2 | ((NUC_U8_TO_BIT_RC[nuc as usize - 65] as u64) << (K-1)*2))
    }
}

#[derive(Clone)]
pub struct KmerIter<'a, const K: usize, const SKIP_N: bool> {
    seq: &'a [u8],
    fwd: Kmer<K>,
    rev: Kmer<K>,
    pos: usize,
    i: usize,
    shift: usize,
    mask: u64,
    init: bool,
    assume_perfect_data: bool,
}

impl<'a, const K: usize, const SKIP_N: bool> KmerIter<'a, K, SKIP_N> {
    pub fn new(seq: &[u8]) -> KmerIter::<K, SKIP_N> {
        KmerIter {
            seq: seq,
            fwd: Kmer::<K>(0),
            rev: Kmer::<K>(0),
            pos: 0,
            i: 0,
            shift: 2 * (K - 1),
            mask: (1 << (2 * K)) - 1,
            init: false,
            assume_perfect_data: true,
        }
    }
    
    pub fn skip_unsupported_nucleotides(&mut self) -> bool {
        assert!(self.pos + K <= self.seq.len());
        // You want to find if there is a nuc not ACGT in the current Kmer slice
        // If pos_skip is called, one of two things happend:
        // 1: It's the first k-mer
        // 2: While constructing the next k-mer, a nuc not ACGT has been encountered
        while let Some(pos) = next_imperfect_nuc(&self.seq[self.pos..self.pos+K]) {
            self.pos = self.pos + pos + 1;
            if self.pos+K >= self.seq.len() {
                return true
            }
        }
        false
    }

    pub fn set(&mut self, seq: &'a[u8]) {
        self.seq = seq
    }

    pub fn set_assume_perfect_data(&mut self, val: bool) {
        self.assume_perfect_data = val;
    }

    pub fn init(&mut self) {
        //let kmer = Kmer::<K>::from_slice(&self.seq[self.pos..self.pos+K]).unwrap();
        self.rev = Kmer::<K>(0);//kmer.rc();
        self.fwd = Kmer::<K>(0);
        self.init = true;
    }

    #[inline]
    pub fn perfect_next_pair(&mut self) -> Option<(usize, Kmer::<K>, Kmer::<K>)> {
        loop {
            let mut nuc: u8 = self.seq[self.pos];
            let allowed = allowed(nuc);
            self.pos += 1;

            if !allowed && !SKIP_N {
                nuc = b'A';
            } 

            match allowed || !SKIP_N {
                false => {
                    self.i = 0;
                    self.fwd.0 = 0;
                    self.rev.0 = 0;
                }
                true => {
                    self.fwd.0 = (self.fwd.0 << 2) | (NUC_U8_TO_BIT[nuc as usize - 65] as u64);
                    self.rev.0 = (self.rev.0 >> 2) | ((3 - NUC_U8_TO_BIT[nuc as usize - 65] as u64) << self.shift);
                    self.i += 1;
                    if self.i == K {
                        self.fwd.0 &= self.mask;
                        self.i -= 1;
                        return Some((self.pos - K, self.fwd.clone(), self.rev.clone()));
                    }
                    if self.pos >= self.seq.len() {
                        return None;
                    }
                }
            }
        }
        None
    }
}


impl<'a, const K: usize, const SKIP_N: bool> Iterator for KmerIter<'a, K, SKIP_N> {
    type Item = (usize, Kmer::<K>, Kmer::<K>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == usize::MAX || self.pos >= self.seq.len() { return None; }
        self.perfect_next_pair()
        // if self.assume_perfect_data {  } else { self.imperfect_next() }
    }
}



#[cfg(test)]
mod kmer_tests {
    use std::{cmp::min, hash::{DefaultHasher, Hash, Hasher}};

    use kmers::KmerIterator;
    use test::Bencher;

    use crate::{
        consecutive::kmer::{Kmer, KmerIter}, 
        kmer_utils::{canonical_rc, canonical_rc2}, 
        minimizer::{context_free::Minimizer, universal_minimizer::{self, UniMin}}, 
        syncmer::closed_syncmer::ClosedSyncmer};
    use kmers::{Kmer as ExtKmer};
    use fxhash::FxHasher64;

    #[test]
    fn test_kmer_iter() {
        const K: usize = 5;
        const SKIP_N: bool = true;
        let kmer: &str = "GACTACTACT";

        let true_kmers = [
            Kmer::<K>::from_slice("GACTA".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("ACTAC".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("CTACT".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("TACTA".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("ACTAC".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("CTACT".as_bytes()).unwrap().get_smallest_rc(),]; 

        let mut count = 0;
        for (pos, fwd, rev) in KmerIter::<K, SKIP_N>::new(kmer.as_bytes()) {
            let kmer = min(fwd, rev);
            assert_eq!(true_kmers[pos].to_string().unwrap(), kmer.to_string().unwrap());
            count += 1;
        }
        assert_eq!(true_kmers.len(), count);
    }

    #[test]
    fn test_kmer_iter_imperfect() {
        const K: usize = 5;
        const SKIP_N: bool = true;
        let kmer: &str = "GANTACTACT";

        let true_kmers = [
            Kmer::<K>::from_slice("TACTA".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("ACTAC".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("CTACT".as_bytes()).unwrap().get_smallest_rc(),]; 

        let mut iter = KmerIter::<K, SKIP_N>::new(kmer.as_bytes());
        iter.set_assume_perfect_data(false);

        let mut count = 0;
        for (_pos, fwd, rev) in iter {
            let kmer = min(fwd, rev);
            assert_eq!(true_kmers[count].to_string().unwrap(), kmer.to_string().unwrap());
            count += 1;
        }
        assert_eq!(true_kmers.len(), count);
    }


    #[test]
    fn test_kmer_iter_imperfect_2() {
        const K: usize = 5;
        let kmer: &str = "GANTACNACT";

        let mut iter = KmerIter::<K, true>::new(kmer.as_bytes());
        iter.set_assume_perfect_data(false);

        let mut count = 0;
        for (_pos, fwd, rev) in iter {
            count += 1;
        }

        assert_eq!(0, count);
    }

    #[test]
    fn test_kmer_iter_imperfect_3() {
        const K: usize = 5;
        let kmer: &str = "GACTACTNCT";

        let true_kmers = [
            Kmer::<K>::from_slice("GACTA".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("ACTAC".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("CTACT".as_bytes()).unwrap().get_smallest_rc(),]; 

        let mut iter = KmerIter::<K, true>::new(kmer.as_bytes());
        iter.set_assume_perfect_data(false);

        let mut count = 0;
        for (_pos, fwd, rev) in iter {
            let kmer = min(fwd, rev);
            assert_eq!(true_kmers[count].to_string().unwrap(), kmer.to_string().unwrap());
            count += 1;
        }

        assert_eq!(3, count);
    }

    #[test]
    fn test_kmer_iter_imperfect_4() {
        const K: usize = 5;
        let kmer: &str = "GACTACTNCTCGTNTCTTC";

        let true_kmers = [
            Kmer::<K>::from_slice("GACTA".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("ACTAC".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("CTACT".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("CTCGT".as_bytes()).unwrap().get_smallest_rc(),
            Kmer::<K>::from_slice("TCTTC".as_bytes()).unwrap().get_smallest_rc(),]; 

        let mut iter = KmerIter::<K, true>::new(kmer.as_bytes());
        iter.set_assume_perfect_data(false);

        let mut count = 0;
        for (_pos, fwd, rev) in iter {
            let kmer = min(fwd, rev);
            assert_eq!(true_kmers[count].to_string().unwrap(), kmer.to_string().unwrap());
            count += 1;
        }

        assert_eq!(5, count);
    }

    #[test]
    fn test_kmer_iter_pos_skip() {
        const K: usize = 5;
        let kmer: &str = "GANTACTACT";

        let mut iter = KmerIter::<K, true>::new(kmer.as_bytes());
        iter.set_assume_perfect_data(false);

        iter.skip_unsupported_nucleotides();
        assert_eq!(iter.pos, 3);
    }

    #[test]
    fn test_kmer_iter_pos_skip2() {
        const K: usize = 5;
        //                0123456789
        let kmer: &str = "GANTACTNCT";

        let mut iter = KmerIter::<K, true>::new(kmer.as_bytes());
        iter.set_assume_perfect_data(false);

        iter.skip_unsupported_nucleotides();
        assert_eq!(iter.pos, 8);
    }


    #[test]
    fn test_create_kmer() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        assert_eq!(
            Kmer::<17>::from_slice(&kmer.as_bytes()).unwrap().0,
            0b1000011100011100011100100111100100
        );
    }

    #[test]
    fn test_create_kmer_2() {
        let kmer: &str = "GACTACTACTAGCTGCA";

        // G A C T A C T A C T A G C T G C A
        // 100001110001110001110010011110
        
        let init: Kmer<15> = Kmer::<15>::from_slice(&kmer[0..15].as_bytes()).unwrap();
        let first = init.append(kmer.as_bytes()[15]);
        let second = first.append(kmer.as_bytes()[16]);


        assert_eq!(
            init.0,
            0b100001110001110001110010011110
        );

        assert_eq!(
            first.0,
            0b000111000111000111001001111001
        );

        assert_eq!(
            second.0,
            0b011100011100011100100111100100
        );
    }


    #[test]
    fn test_cut_middle() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        assert_eq!(
            Kmer::<17>::from_slice(&kmer.as_bytes()).unwrap().0,
            0b1000011100011100011100100111100100
        );
    }

    //----------------------- BENCHMARKS ------------------------------------

    #[bench]
    fn bench_construct_kmer(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCA";
        b.iter(|| Kmer::<17>::from_slice(&kmer.as_bytes()))
    }

    #[bench]
    fn bench_construct_kmer_external(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCA";
        b.iter(|| ExtKmer::make(kmer));
    }

    #[bench]
    fn bench_construct_kmers(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGC";
        const K: usize = 15; 
        let mut sum = 0;
        b.iter(|| {
            let init: Kmer<K> = Kmer::<K>::from_slice(&kmer[0..K].as_bytes()).unwrap();
            sum += init.0;
            for i in K..kmer.len() {
                let first = init.append(kmer.as_bytes()[i]);
                sum += first.0;
            }
        })
    }

    #[bench]
    fn bench_construct_kmers_naive(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        let mut sum = 0;
        b.iter(|| {
            for (start, end) in (K..kmer.len()).enumerate() {
                let next = Kmer::<K>::from_slice(&kmer[start..end].as_bytes()).unwrap();
                sum += min(next.0, canonical_rc2::<K>(next.0));
            }
        })
    }

    #[bench]
    fn bench_construct_kmers_iterator(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        let mut sum = 0;
        b.iter(|| {
            KmerIter::<K, true>::new(&kmer.as_bytes()).for_each(|(_, fwd, rev)| {
                let kmer = min(fwd, rev);
                sum += kmer.0;
            });
        })
    }

    #[bench]
    fn bench_construct_kmers_iterator_option(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        let mut sum = 0;
        b.iter(|| {
            let mut k_iter = KmerIter::<K, true>::new(&kmer.as_bytes());
            k_iter.set_assume_perfect_data(false);
            k_iter.for_each(|(_, fwd, rev)| {
                let kmer = min(fwd, rev);
                sum += kmer.0;
            });
        })
    }

    #[bench]
    fn bench_construct_kmers_iterator_option_perfect(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        let mut sum = 0;
        b.iter(|| {
            let mut k_iter = KmerIter::<K, true>::new(&kmer.as_bytes());
            k_iter.set_assume_perfect_data(true);
            for (pos, fwd, rev) in k_iter {
                let kmer = min(fwd, rev);
                sum += kmer.0;
            }
        })
    }

    #[bench]
    fn bench_construct_kmers_iterator_option_closed_syncmer(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        const S: usize = 7; 
        let mut cs = ClosedSyncmer::<{K as u64},{S as u64}>::new();

        let mut sum = 0;
        b.iter(|| {
            let mut k_iter = KmerIter::<K, true>::new(&kmer.as_bytes());
            k_iter.set_assume_perfect_data(false);

            k_iter.for_each(|(_, fwd, rev)| {
                let kmer = min(fwd, rev);
                sum += kmer.0;
                sum += cs.is_minimizer(kmer.0) as u64;
            });
        })
    }

    #[bench]
    fn bench_construct_kmers_iterator_option_unimin(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        
        const K: usize = 15; 
        const S: usize = 7; 
        let mut minim: UniMin<15, 8> = UniMin::<K,8>::new();

        // let mut state = DefaultHasher::new();
        let mut hasher = FxHasher64::default();

        let mut kmers: Vec<(usize, Kmer<15>)> = vec![(0, Kmer::<K>(0)); 50];

        let mut sum = 0;
        b.iter(|| {
            let mut k_iter = KmerIter::<K, true>::new(&kmer.as_bytes());

            k_iter.for_each(|(_, fwd, rev)| {
                let kmer = min(fwd, rev);
                kmer.hash(&mut hasher);
                let hash = hasher.finish();
                sum += minim.is_minimizer(hash >> (64-(K*2))) as usize;
            });
        })
    }

    #[bench]
    fn bench_construct_kmers_external_lib(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        

        const K: usize = 15; 
        const S: usize = 7; 
        let mut cs = ClosedSyncmer::<{K as u64},{S as u64}>::new();

        let mut sum = 0;
        b.iter(|| {
            let kit = KmerIterator::new(K, kmer.as_bytes().iter());
            for (kmer1, kmer2) in kit {
                sum += min(kmer1.0, kmer2.0);
            }
        });
    }
    #[bench]
    fn bench_construct_kmers_external_lib2(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCGACTACTACTAGCTGCAGTACGTACGTCATCGTCTCTCGTGCTCCCGTGCCCGTGC";
        

        const K: usize = 15;

        let mut sum = 0;

        b.iter(|| {
            ExtKmer::with_many(K, &kmer.as_bytes().iter(), |kmer| {
                sum += kmer.0;
                sum += kmer.rev(K).0;
            });
        });
    }
}