use std::cmp::min;

use crate::{minimizer::context_free::Minimizer, utils::{min_index}};

pub struct ClosedSyncmer<const K: u64, const S: u64>
        where [(); (K - S + 1) as usize]: {
    smers_count: usize,
    smers: [u64; (K - S + 1) as usize] ,
    shifts: [u64; (K - S + 1) as usize] ,
    bitmask: u64,
}

impl<const K: u64, const S: u64> ClosedSyncmer<K, S> where [(); (K - S + 1) as usize]: {
    pub fn new() -> Self {
        let smers_count = (K - S + 1) as usize;
        let mut ret_obj = Self {
            smers_count: smers_count,
            smers: [0; (K - S + 1) as usize],
            shifts: [0; (K - S + 1) as usize],
            bitmask: (1 << 2*S) - 1,
        };
        for i in 0..ret_obj.smers_count {
            ret_obj.shifts[i] = ((smers_count*2) - ((i+1)*2)) as u64;
        }

        ret_obj
    }

    fn load(&mut self, canonical_kmer: u64) {
        for i in 0..self.smers_count {
            self.smers[i] = (canonical_kmer >> self.shifts[i]) & self.bitmask;
        }
    }

    pub fn index_min(&self) -> usize {
        return min_index(self.smers).unwrap();
    }
}

impl<const K: u64, const S: u64> Minimizer for ClosedSyncmer <K, S> where [(); (K - S + 1) as usize]: {
    fn is_minimizer(&mut self, hash: u64) -> bool {
        if S == 0 {
            return false;
        }
        self.load(hash);
        let min = unsafe { self.smers.iter().min_by(|a, b| a.partial_cmp(b).unwrap_unchecked())}.unwrap();
        min == &self.smers[0] || min == self.smers.last().unwrap()
    }
}


#[cfg(test)]
mod closed_syncmer_tests {
    use std::cmp::min;

    use test::Bencher;

    use crate::consecutive::kmer::{Kmer, KmerIter};
    use crate::minimizer::context_free::Minimizer;
    use crate::syncmer::closed_syncmer::ClosedSyncmer;
    use crate::kmer_utils::{self, canonical_rc, slice_to_canonical_4};

    #[test]
    fn test_closed_syncmer_load_1() {
        let mut cs = ClosedSyncmer::<15,7>::new();

        let kmer = "ACGTGTCTCGATAGC";

        let fwd = slice_to_canonical_4(kmer.as_bytes());
        assert_eq!(fwd, min(fwd, canonical_rc::<15>(fwd)));

        assert_eq!(cs.bitmask, 0b11111111111111);

        cs.load(fwd);
        //  A C G T G T C T C G A T A G C
        // 00011011101101
        //   01101110110111
        //     10111011011101
        //       11101101110110
        //         10110111011000
        //           11011101100011
        //             01110110001100
        //               11011000110010
        //                 01100011001001
        let correct = [
            0b00011011101101,
            0b01101110110111,
            0b10111011011101,
            0b11101101110110,
            0b10110111011000,
            0b11011101100011,
            0b01110110001100,
            0b11011000110010,
            0b01100011001001,
        ];
        for i in 0..7 {
            assert_eq!(cs.smers[i], correct[i]);
        }
        assert_eq!(cs.index_min(), 0);
        assert!(cs.is_minimizer(fwd));
    }
}