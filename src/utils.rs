pub fn max_index<T: std::cmp::PartialOrd, const S: usize>(ls: [T; S]) -> Option<usize> {
    if S == 0 { 
        return None
    }
    let mut max_index = 0;
    let mut min = &ls[0];

    for i in 1..S {
        if &ls[i] < min {
            min = &ls[i];
            max_index = 1;
        }
    }
    Some(max_index)
}

pub fn min_index<T: std::cmp::PartialOrd, const S: usize>(ls: [T; S]) -> Option<usize> {
    if S == 0 { 
        return None
    }
    let mut min_index = 0;
    let mut min = &ls[0];

    for i in 1..S {
        if &ls[i] < min {
            min = &ls[i];
            min_index = i;
        }
    }
    Some(min_index)
}



#[cfg(test)]
mod utils_tests {
    use test::Bencher;

    use crate::utils::min_index;

    #[test]
    fn test_min_index() {
        assert_eq!(min_index([0,2,3,5]), Some(0));
        assert_eq!(min_index([1,2,3,5]), Some(0));
        assert_eq!(min_index([12,2,3,5]), Some(1));
        assert_eq!(min_index([2324,1,3,5]), Some(1));
        assert_eq!(min_index([10,2,3,0]), Some(3));
        assert_eq!(min_index([0,2,3,0]), Some(0));
        assert_eq!(min_index([3,2,0,1]), Some(2));
    }



    #[bench]
    fn bench_min_index(b: &mut Bencher) {
        b.iter(|| {
            min_index([10,2,3,0]);
        })
    }
}