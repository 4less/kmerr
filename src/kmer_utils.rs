

/// Fastest way to convert u8 into 2-bit canonical representation
/// 
pub static NUC_U8_TO_BIT: [u8; 20] = [
    0b00, 0b111, 0b01, 0b111, 0b111, 0b111, 0b10, 0b111, 0b111, 0b111, 0b111, 0b111, 0b111, 0b111,
    0b111, 0b111, 0b111, 0b111, 0b111, 0b11,
];
pub static NUC_U8_TO_BIT_RC: [u8; 20] = [
    0b11, 0b111, 0b10, 0b111, 0b111, 0b111, 0b01, 0b111, 0b111, 0b111, 0b111, 0b111, 0b111, 0b111,
    0b111, 0b111, 0b111, 0b111, 0b111, 0b00,
];

pub static NUC_BIT_RC: [u8; 4] = [
    0b11, 0b10, 0b01, 0b00,
];

pub static NUC_TABLE_REV: [u8; 4] = [
    b'A', b'C', b'G', b'T',
];

pub fn u8_to_canonical(c: u8) -> u8 {
    match c {
        b'A' => 0b00,
        b'C' => 0b01,
        b'G' => 0b10,
        b'T' => 0b11,
        _ => 0b00,
    }
}

pub fn slice_to_canonical(slice: &[u8]) -> u64 {
    if slice.len() > 32 {
        panic!("u64 can only hold 32 bases")
    };
    let mut canon: u64 = 0;
    for (i, &c) in slice.iter().rev().enumerate() {
        canon |= (u8_to_canonical(c) as u64) << i * 2;
    }
    canon
}

pub fn canonical_to_u8<const K: usize>(mut canon: u64) -> [u8; K] {
    let mut res= [0u8; K];
    for i in (0..K).rev() {
        res[i] = NUC_TABLE_REV[canon as usize & 0b11];
        canon >>= 2;
    }
    res
}

pub fn nuc_perfect(c: u8) -> bool {
    c == b'A' || c == b'C' || c == b'G' || c == b'T'
}

pub fn next_imperfect_nuc(seq: &[u8]) -> Option<usize> {
    seq.iter().position(|&c| {
        !nuc_perfect(c)
    })
}

pub fn canonical_rc<const K: usize>(mut canon: u64) -> u64 {
    let mut canon_rc = NUC_BIT_RC[(canon & 0b11) as usize] as u64;
    for _ in 1..K {
        canon_rc <<= 2;
        canon >>= 2;
        canon_rc |= NUC_BIT_RC[(canon & 0b11) as usize] as u64;
    }
    canon_rc
}

pub fn canonical_rc2<const K: usize>(canon: u64) -> u64 {
    rev::<K>(!canon)
}

/// Create a reversed k-mer.
#[inline]
pub fn rev<const K: usize>(mut x: u64) -> u64 {
    const M2: u64 = 0x3333333333333333;
    const M3: u64 = 0x0F0F0F0F0F0F0F0F;
    const M4: u64 = 0x00FF00FF00FF00FF;
    const M5: u64 = 0x0000FFFF0000FFFF;
    const M6: u64 = 0x00000000FFFFFFFF;

    x = ((x >> 2) & M2) | ((x & M2) << 2);
    x = ((x >> 4) & M3) | ((x & M3) << 4);
    x = ((x >> 8) & M4) | ((x & M4) << 8);
    x = ((x >> 16) & M5) | ((x & M5) << 16);
    x = ((x >> 32) & M6) | ((x & M6) << 32);
    x >>= 64 - 2 * K;
    x
}

pub fn slice_to_canonical_2(slice: &[u8]) -> u64 {
    if slice.len() > 32 {
        panic!("u64 can only hold 32 bases")
    };
    slice.iter().rev().enumerate().fold(0u64, |acc, (i, &c)| {
        acc | (u8_to_canonical(c) as u64) << i * 2
    })
}

pub fn slice_to_canonical_3(slice: &[u8]) -> u64 {
    if slice.len() > 32 {
        panic!("u64 can only hold 32 bases")
    };
    slice.iter().rev().enumerate().fold(0u64, |acc, (i, &c)| {
        acc | (NUC_U8_TO_BIT[c as usize - 65] as u64) << i * 2
    })
}

/// Fastest, see benchmarks
pub fn slice_to_canonical_4(slice: &[u8]) -> u64 {
    if slice.len() > 32 {
        panic!("u64 can only hold 32 bases")
    };
    let mut acc = 0;
    let mut i  = 0;
    for &c in slice {
        i += 1;
        acc |= (NUC_U8_TO_BIT[c as usize - 65] as u64) << (slice.len() - i) * 2;
    }
    acc
}

#[inline]
pub const fn allowed(nuc: u8) -> bool {
    nuc == b'A' || nuc == b'C' || nuc == b'G' || nuc == b'T'
}

pub fn find_init_position(seq: &[u8]) -> Option<usize> {
    let mut counter = 0;
    let mut counter_start = 0;
    for (i, &c) in seq.iter().enumerate() {
        if counter == 15 { return Some(counter_start) }
        counter += 1;
        if !allowed(c) {
            counter = 0;
            counter_start = i+1;
        }
    };
    None
}

pub fn hamming_distance(_a: u64, _b: u64) {
    // (a & b).popcnt()
}


#[cfg(test)]
mod kmer_utils_tests {
    use test::Bencher;

    use crate::kmer_utils::{canonical_rc, canonical_rc2, canonical_to_u8, hamming_distance, slice_to_canonical, slice_to_canonical_3, slice_to_canonical_4, u8_to_canonical};

    use super::slice_to_canonical_2;

    #[test]
    fn test_char_to_canonical1() {
        assert_eq!(u8_to_canonical(b'A'), 0b00);
    }

    #[test]
    fn test_char_to_canonical2() {
        assert_eq!(u8_to_canonical(b'Z'), 0b00);
    }

    #[test]
    fn test_char_to_canonical3() {
        assert_eq!(u8_to_canonical(b'T'), 0b11);
    }

    #[test]
    fn test_slice_to_canonical1_1() {
        let kmer: &str = "ACGT";
        assert_eq!(slice_to_canonical(&kmer.as_bytes()), 0b00011011);
    }

    #[test]
    fn test_slice_to_canonical1_2() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        assert_eq!(
            slice_to_canonical(&kmer.as_bytes()),
            0b1000011100011100011100100111100100
        );
    }

    #[test]
    fn test_slice_to_canonical2_1() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        assert_eq!(
            slice_to_canonical_2(&kmer.as_bytes()),
            0b1000011100011100011100100111100100
        );
    }

    #[test]
    fn test_slice_to_canonical3() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        assert_eq!(
            slice_to_canonical_3(&kmer.as_bytes()),
            0b1000011100011100011100100111100100
        );
    }

    #[test]
    fn test_slice_to_canonical_4() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        assert_eq!(
            slice_to_canonical_4(&kmer.as_bytes()),
            0b1000011100011100011100100111100100
        );
    }

    #[test]
    fn test_canonical_to_u8_1() {
        let kmer: &str = "GACTACTACTAGCTGCA";
        let canon = canonical_to_u8::<17>(slice_to_canonical_4(&kmer.as_bytes()));
        assert_eq!(
            std::str::from_utf8(&canon).unwrap(),
            kmer
        );
    }


    #[test]
    fn test_canonical_rc_1() {
        let kmer_fwd: &str = "GACTACTACTAGCCA";
        let kmer_rev: &str = "TGGCTAGTAGTAGTC";
        let fwd = slice_to_canonical_4(kmer_fwd.as_bytes());
        let rev: u64 = slice_to_canonical_4(kmer_rev.as_bytes());
        assert_eq!(fwd, canonical_rc::<15>(rev));
        assert_eq!(rev, canonical_rc::<15>(fwd));
        assert_ne!(fwd, rev);
    }

    #[test]
    fn test_canonical_rc_2() {
        let kmer_fwd: &str = "GACTACTACTAGCCA";
        let kmer_rev: &str = "TGGCTAGTAGTAGTC";
        let fwd = slice_to_canonical_4(kmer_fwd.as_bytes());
        let rev: u64 = slice_to_canonical_4(kmer_rev.as_bytes());
        assert_eq!(fwd, canonical_rc2::<15>(rev));
        assert_eq!(rev, canonical_rc2::<15>(fwd));
        assert_ne!(fwd, rev);
    }

    // #[test]
    // fn test_slice_to_canonical_4() {
    //     let a = 0b0100011000101110101011010101000;
    //     let b = 0b0101011000101010001011110101110;
    //     assert_eq!(
    //         hamming_distance(a, b),
    //         6
    //     );
    // }


    //----------------------- BENCHMARKS ------------------------------------

    #[bench]
    fn bench_slice_to_canonical1_1(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCA";
        b.iter(|| slice_to_canonical(kmer.as_bytes()))
    }

    #[bench]
    fn bench_slice_to_canonical2_1(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCA";
        b.iter(|| slice_to_canonical_2(kmer.as_bytes()))
    }

    #[bench]
    fn bench_slice_to_canonical3_1(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCA";
        b.iter(|| slice_to_canonical_3(kmer.as_bytes()))
    }

    #[bench]
    fn bench_slice_to_canonical4_1(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCTGCA";
        b.iter(|| slice_to_canonical_4(kmer.as_bytes()))
    }

    #[bench]
    fn bench_canonical_rc_1(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCCA";
        let fwd = slice_to_canonical_4(kmer.as_bytes());
        b.iter(|| canonical_rc::<15>(fwd))
    }


    #[bench]
    fn bench_canonical_rc_2(b: &mut Bencher) {
        let kmer: &str = "GACTACTACTAGCCA";
        let fwd = slice_to_canonical_4(kmer.as_bytes());
        b.iter(|| canonical_rc2::<15>(fwd))
    }
}
