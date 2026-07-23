//! Uniformity check for core-mer key placement under different masking schemes.
//!
//! Sweeps the entire `2^(2C)` core-mer space, applies canonicalisation + the closed-syncmer
//! selection, and histograms where the *selected* core-mers' keys land across the range under four
//! schemes. The key distribution decides how balanced the physical shard slices can be.
//!
//! Schemes:
//!   * `before`       -- current: syncmer(canonical), key = canonical.
//!   * `mask_syncmer` -- mask before the syncmer, key = mask(canonical). (mask feeds the selection)
//!   * `mask_key`     -- syncmer(canonical), key = mask(canonical).       (mask only the key)
//!   * `mask_both`    -- syncmer(maskA(canonical)), key = maskB(canonical). (two independent hashes)
//!
//! The point it demonstrates: the closed-syncmer rule keeps a value iff its *leading* s-mer is the
//! minimum, which biases the kept *value* toward small numbers. So a key that is also the syncmer's
//! input (`before`, `mask_syncmer`) stays piled up at the low end. Only when the key is a hash that
//! is independent of the syncmer's input (`mask_key`, `mask_both`) does it spread uniformly.
//!
//! `#[ignore]`d because it visits ~10^9 core-mers; run it in release:
//!
//! ```text
//! cargo +nightly test --release -p kmerrs mask_uniformity -- --ignored --nocapture
//! # optional: COREMER_DENSITY_OUT=/path/density.tsv to choose the TSV location
//! ```

use std::io::Write;

use crate::kmer_utils::{canonical_rc, mask_coremer};
use crate::minimizer::context_free::Minimizer;
use crate::syncmer::closed_syncmer::ClosedSyncmer;

const C: usize = 15; // core-mer length (matches flexalign)
const S: usize = 4; // s-mer length for the closed syncmer
const L: usize = C - S + 1; // number of s-mers per core-mer
const BITS: u32 = 2 * C as u32; // 30
const KEY_SPACE: u64 = 1 << BITS; // 2^30
const LOG_BINS: u32 = 10;
const NBINS: usize = 1 << LOG_BINS; // 1024 bins over the range
const BIN_SHIFT: u32 = BITS - LOG_BINS; // key >> 20 -> [0, 1024)

/// A second bijection over the `2*C`-bit space, independent of [`mask_coremer`] (different odd
/// multipliers / shifts). Used as the key hash in the two-hash scheme so the key does not correlate
/// with the syncmer's masked input.
#[inline(always)]
fn mask_coremer_b<const C: usize>(x: u64) -> u64 {
    let mask = (1u64 << (2 * C)) - 1;
    let mut h = x & mask;
    h ^= h >> 16;
    h = h.wrapping_mul(0x27D4EB2F) & mask;
    h ^= h >> 11;
    h = h.wrapping_mul(0x165667B1) & mask;
    h ^= h >> 15;
    h & mask
}

/// Density divisor for the modulo sampler (keep iff `mask(coremer) % MOD == 0`); ~1/MOD of values.
const MOD: u64 = 10;

const NS: usize = 7;
const SCHEMES: [&str; NS] = [
    "syncmer(uniform)", // pure test: is_minimizer over every 30-bit value (perfect uniform input)
    "canonical(all)",
    "before",
    "mask_syncmer",
    "mask_key",
    "mask_both",
    "mask_modulo", // single mask + uniform-density sampler: sel = mask%MOD==0, key = mask
];

struct Hist {
    bins: [Vec<u64>; NS], // one histogram per SCHEMES entry
}

impl Hist {
    fn zero() -> Self {
        Self { bins: std::array::from_fn(|_| vec![0u64; NBINS]) }
    }
    fn merge(&mut self, o: &Hist) {
        for s in 0..NS {
            for i in 0..NBINS {
                self.bins[s][i] += o.bins[s][i];
            }
        }
    }
}

/// Sweep `[lo, hi)`, counting each canonical core-mer exactly once (only when `x` is the canonical
/// representative of its {x, rc(x)} pair).
fn sweep(lo: u64, hi: u64) -> Hist {
    let mut cs = ClosedSyncmer::<C, S, L>::new();
    let mut h = Hist::zero();
    let bin = |k: u64| (k >> BIN_SHIFT) as usize;
    for x in lo..hi {
        // Scheme 0: feed the syncmer a perfectly uniform input (every value once) and record which
        // values it keeps. This isolates the selection rule from any masking/canonicalisation.
        if cs.is_minimizer(x) {
            h.bins[0][bin(x)] += 1;
        }

        let rc = canonical_rc::<C>(x);
        if x > rc {
            continue; // rc is the canonical rep; counted when we reach it
        }
        let canonical = x; // == min(x, rc)
        let ma = mask_coremer::<C>(canonical);
        let mb = mask_coremer_b::<C>(canonical);
        let sel_plain = cs.is_minimizer(canonical);
        let sel_masked = cs.is_minimizer(ma);

        h.bins[1][bin(canonical)] += 1; // all canonical (no selection)
        if sel_plain {
            h.bins[2][bin(canonical)] += 1; // before:       key = canonical
            h.bins[4][bin(ma)] += 1; //         mask_key:     key = mask(canonical)
        }
        if sel_masked {
            h.bins[3][bin(ma)] += 1; //         mask_syncmer: key = mask(canonical)
            h.bins[5][bin(mb)] += 1; //         mask_both:    key = maskB(canonical)
        }
        if ma % MOD == 0 {
            h.bins[6][bin(ma)] += 1; //         mask_modulo:  key = mask(canonical)
        }
    }
    h
}

fn parallel_histogram() -> Hist {
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8).max(1);
    let chunk = KEY_SPACE / threads as u64;
    std::thread::scope(|sc| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let lo = t as u64 * chunk;
                let hi = if t == threads - 1 { KEY_SPACE } else { lo + chunk };
                sc.spawn(move || sweep(lo, hi))
            })
            .collect();
        let mut total = Hist::zero();
        for handle in handles {
            total.merge(&handle.join().unwrap());
        }
        total
    })
}

/// sum, min, max, and coefficient of variation (std/mean) over a histogram.
fn stats(h: &[u64]) -> (u64, u64, u64, f64) {
    let sum: u64 = h.iter().sum();
    let n = h.len() as f64;
    let mean = sum as f64 / n;
    let min = *h.iter().min().unwrap();
    let max = *h.iter().max().unwrap();
    let var = h.iter().map(|&c| (c as f64 - mean).powi(2)).sum::<f64>() / n;
    (sum, min, max, var.sqrt() / mean.max(1.0))
}

/// Compact ASCII bar chart: fold the NBINS bins down to `cols` columns of block glyphs.
fn ascii_plot(label: &str, h: &[u64], cols: usize) {
    let group = h.len() / cols;
    let folded: Vec<u64> = (0..cols).map(|c| h[c * group..(c + 1) * group].iter().sum()).collect();
    let peak = *folded.iter().max().unwrap().max(&1);
    let glyphs = [' ', '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    let bar: String = folded
        .iter()
        .map(|&v| glyphs[((v as f64 / peak as f64) * 8.0).round() as usize])
        .collect();
    println!("{label:<15} |{bar}|");
}

#[test]
#[ignore = "sweeps all 2^30 core-mers; run with --release --ignored --nocapture"]
fn mask_uniformity() {
    eprintln!("sweeping {KEY_SPACE} core-mers ({NBINS} bins) ...");
    let t0 = std::time::Instant::now();
    let h = parallel_histogram();
    eprintln!("swept in {:?}", t0.elapsed());

    let all_total: u64 = h.bins[1].iter().sum(); // canonical(all)
    println!("\n=== core-mer key density across [0, 2^{BITS}) ({NBINS} bins) ===");
    for s in 0..NS {
        ascii_plot(SCHEMES[s], &h.bins[s], 96);
    }
    println!("\n{:<16} {:>12} {:>8} {:>10} {:>10} {:>8}", "scheme", "selected", "density", "bin_min", "bin_max", "cv");
    for s in 0..NS {
        let (sum, min, max, cv) = stats(&h.bins[s]);
        println!(
            "{:<16} {sum:>12} {:>8.4} {min:>10} {max:>10} {cv:>8.4}",
            SCHEMES[s],
            sum as f64 / all_total as f64
        );
    }

    // Write the full histogram as TSV for external plotting.
    let out = std::env::var("COREMER_DENSITY_OUT").unwrap_or_else(|_| "coremer_density.tsv".into());
    let mut f = std::fs::File::create(&out).expect("create density tsv");
    writeln!(f, "bin\tkey_lo\t{}", SCHEMES.join("\t")).unwrap();
    for i in 0..NBINS {
        write!(f, "{i}\t{}", (i as u64) << BIN_SHIFT).unwrap();
        for s in 0..NS {
            write!(f, "\t{}", h.bins[s][i]).unwrap();
        }
        writeln!(f).unwrap();
    }
    println!("\nwrote {out}");

    // The schemes whose key is independent of the selection statistic (mask_key=4, mask_both=5,
    // mask_modulo=6) must be ~uniform: no empty bins and a small coefficient of variation.
    for &s in &[4usize, 5, 6] {
        let (_, min, _, cv) = stats(&h.bins[s]);
        assert!(min > 0, "{}: has an empty bin", SCHEMES[s]);
        assert!(cv < 0.05, "{}: not uniform (cv {cv:.4} >= 0.05)", SCHEMES[s]);
    }
}
