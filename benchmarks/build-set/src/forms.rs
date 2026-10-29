//! The draw behind a `forms` set: a source sampled to equal counts per
//! doubling length bin, and written three ways, as drawn, backwards, and with
//! each sequence's residues shuffled.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, ensure};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use libsail::collection::{Aggregate, Indexable};
use libsail::select::Selection;
use libsail::seq::fasta::{self, DEFAULT_LINE_WIDTH, IndexedFasta};

/// The three ways the draw is written, in the order they are written.
pub const FORMS: [&str; 3] = ["original", "reversed", "shuffled"];

/// The shortest sequence drawn, and the lower bound of the first bin.
pub const SHORTEST: usize = 16;

/// The number of length bins: `[16, 32)` doubling up to an open `[16384, inf)`.
pub const BINS: usize = 11;

/// What `bins` records for a sequence shorter than [`SHORTEST`].
const UNDRAWN: u8 = u8::MAX;

/// Records read before the shuffling and wrapping go wide.
const BATCH: usize = 16384;

/// The bin a sequence of `len` residues falls in, or `None` under [`SHORTEST`].
pub fn bin_of(len: usize) -> Option<usize> {
    if len < SHORTEST {
        return None;
    }

    // 16..=31 -> 1 -> 0, 32..=63 -> 2 -> 1, and so on
    Some(((len / SHORTEST).ilog2() as usize).min(BINS - 1))
}

/// The shortest length bin `bin` holds.
pub fn lower(bin: usize) -> usize {
    SHORTEST << bin
}

/// The bin of every record in `seqs`, in collection order, with [`UNDRAWN`]
/// for the ones too short to draw.
pub fn bins(seqs: &Aggregate<IndexedFasta>) -> anyhow::Result<Vec<u8>> {
    let mut out = Vec::with_capacity(seqs.len());

    for part in seqs.parts() {
        // the residues are counted off the framed bytes rather than a parsed
        // record, which would allocate a name and a sequence per record
        let lens = (0..part.len())
            .into_par_iter()
            .map(|n| {
                let bytes = part
                    .record(n)?
                    .with_context(|| format!("record {n} is past the index"))?;

                Ok(bin_of(fasta::len_of(&bytes)).map_or(UNDRAWN, |b| b as u8))
            })
            .collect::<anyhow::Result<Vec<u8>>>()?;

        out.extend(lens);
    }

    Ok(out)
}

/// How many records each bin holds.
pub fn counts(bins: &[u8]) -> [usize; BINS] {
    let mut counts = [0; BINS];
    for &b in bins.iter().filter(|&&b| b != UNDRAWN) {
        counts[b as usize] += 1;
    }
    counts
}

/// Compute each bin's quota: an equal share of `n`, where a bin holding less
/// than its share gives all it has and the rest is split over the others.
///
/// The quotas sum to `n` unless the bins hold fewer than `n` between them, in
/// which case every bin gives all it has.
pub fn quotas(counts: &[usize; BINS], n: usize) -> [usize; BINS] {
    let mut quota = [0; BINS];

    let mut by_size: Vec<usize> = (0..BINS).collect();
    by_size.sort_by_key(|&b| (counts[b], b));

    let mut left = n;
    for (i, &b) in by_size.iter().enumerate() {
        let open = BINS - i;
        let share = left / open;

        if counts[b] <= share {
            quota[b] = counts[b];
            left -= counts[b];
            continue;
        }

        // every bin from here on holds more than the share, since they are in
        // ascending order of size, so each takes the share and the remainder
        // goes one apiece to the shortest of them
        let mut rest: Vec<usize> = by_size[i..].to_vec();
        rest.sort();
        for (j, &r) in rest.iter().enumerate() {
            quota[r] = share + usize::from(j < left % open);
        }
        break;
    }

    quota
}

/// The records to draw: `quotas[b]` of bin `b`'s, uniformly within each bin,
/// as ascending positions in the collection.
pub fn pick(bins: &[u8], quotas: &[usize; BINS], seed: u64) -> Vec<usize> {
    let counts = counts(bins);

    // which members of each bin are kept, counted from the bin's first record;
    // a seed per bin so one bin's draw does not move with another's size
    let mut chosen: Vec<std::iter::Peekable<Selection>> = (0..BINS)
        .map(|b| Selection::new(counts[b], quotas[b], seed.wrapping_add(b as u64)).peekable())
        .collect();
    let mut rank = [0usize; BINS];

    let mut out = Vec::with_capacity(quotas.iter().sum());
    for (n, &b) in bins.iter().enumerate() {
        if b == UNDRAWN {
            continue;
        }
        let b = b as usize;

        if chosen[b].peek() == Some(&rank[b]) {
            chosen[b].next();
            out.push(n);
        }
        rank[b] += 1;
    }

    out
}

/// What one form's file came to. The three forms hold the same lengths under
/// the same names, so one count describes each of them.
pub struct Written {
    pub seqs: usize,
    pub residues: u64,
    pub bytes: u64,
}

/// Write the records at `picked` to `<form>.fa` in `dir` for each of [`FORMS`].
///
/// `bins` is checked against what each record parses to, so a length counted
/// off the framed bytes that disagrees with the parser fails the draw.
pub fn write(
    seqs: &Aggregate<IndexedFasta>,
    bins: &[u8],
    picked: &[usize],
    seed: u64,
    dir: &Path,
) -> anyhow::Result<Written> {
    std::fs::create_dir_all(dir).with_context(|| format!("failed to create {}", dir.display()))?;

    let mut writers = FORMS
        .iter()
        .map(|form| {
            let path = dir.join(format!("{form}.fa"));
            let file = File::create(&path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            Ok(BufWriter::with_capacity(1 << 20, file))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let mut residues = 0u64;

    for batch in picked.chunks(BATCH) {
        let formatted: Vec<(u64, [Vec<u8>; 3])> = batch
            .par_iter()
            .map(|&n| {
                let mut rec = seqs
                    .get(n)
                    .with_context(|| format!("failed to read record {n}"))?;

                ensure!(
                    bin_of(rec.len()).map_or(UNDRAWN, |b| b as u8) == bins[n],
                    "record {n} parses to {} residues, outside the bin it was counted in",
                    rec.len()
                );

                let mut out: [Vec<u8>; 3] = Default::default();
                rec.write_to(&mut out[0], DEFAULT_LINE_WIDTH)?;

                let original = rec.seq.clone();
                rec.reverse();
                rec.write_to(&mut out[1], DEFAULT_LINE_WIDTH)?;

                // seeded by the record's position in the source rather than
                // by its place in the draw, so a sequence shuffles the same
                // way whatever else was drawn beside it
                rec.seq = original;
                let mut rng = StdRng::seed_from_u64(seed ^ (n as u64).wrapping_mul(SPREAD));
                rec.seq.shuffle(&mut rng);
                rec.write_to(&mut out[2], DEFAULT_LINE_WIDTH)?;

                Ok((rec.len() as u64, out))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        for (len, out) in formatted {
            residues += len;
            for (w, bytes) in writers.iter_mut().zip(&out) {
                w.write_all(bytes)?;
            }
        }
    }

    for w in &mut writers {
        w.flush()?;
    }

    let first = dir.join(format!("{}.fa", FORMS[0]));
    Ok(Written {
        seqs: picked.len(),
        residues,
        bytes: std::fs::metadata(&first)
            .with_context(|| format!("failed to stat {}", first.display()))?
            .len(),
    })
}

/// The 64-bit golden ratio, which spreads consecutive positions across the
/// seed space before they are mixed with the recipe's seed.
const SPREAD: u64 = 0x9E37_79B9_7F4A_7C15;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_double_from_sixteen_and_the_last_is_open() {
        assert_eq!(bin_of(15), None);
        assert_eq!(bin_of(16), Some(0));
        assert_eq!(bin_of(31), Some(0));
        assert_eq!(bin_of(32), Some(1));
        assert_eq!(bin_of(16383), Some(9));
        assert_eq!(bin_of(16384), Some(10));
        assert_eq!(bin_of(1_000_000), Some(10));
        assert_eq!(lower(10), 16384);
    }

    #[test]
    fn a_short_bin_hands_its_quota_to_the_rest() {
        let mut counts = [1000; BINS];
        counts[10] = 3;
        counts[9] = 50;

        let q = quotas(&counts, 1100);
        assert_eq!(q.iter().sum::<usize>(), 1100);
        assert_eq!(q[10], 3);
        assert_eq!(q[9], 50);

        // 1047 over nine bins is 116 with 3 over, given to the shortest
        assert_eq!(&q[..9], &[117, 117, 117, 116, 116, 116, 116, 116, 116]);
    }

    #[test]
    fn a_source_smaller_than_the_draw_gives_everything() {
        let counts = [5; BINS];
        assert_eq!(quotas(&counts, 1000), [5; BINS]);
    }

    #[test]
    fn a_pick_meets_every_quota_in_ascending_order() {
        // record n sits in bin n % BINS, and every fifth is too short
        let bins: Vec<u8> = (0..2000)
            .map(|n| match n % 5 {
                0 => UNDRAWN,
                _ => (n % BINS) as u8,
            })
            .collect();

        let q = quotas(&counts(&bins), 500);
        let picked = pick(&bins, &q, 67779);

        assert_eq!(picked.len(), 500);
        assert!(picked.windows(2).all(|w| w[0] < w[1]));

        let mut got = [0; BINS];
        for &n in &picked {
            assert_ne!(bins[n], UNDRAWN);
            got[bins[n] as usize] += 1;
        }
        assert_eq!(got, q);

        assert_eq!(picked, pick(&bins, &q, 67779));
        assert_ne!(picked, pick(&bins, &q, 1));
    }
}
