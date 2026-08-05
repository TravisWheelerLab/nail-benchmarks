//! Fasta shaping used to construct this benchmark. These were the standalone
//! `fab`, `far`, and `fas` binaries; nothing outside mgnify used them.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use bioio::fasta::FastaByteIndex;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;

/// Deal a fasta into `n` shards, shuffling the destination order every `n`
/// records so shards stay comparable in composition rather than reflecting
/// whatever order the source happened to be in.
pub fn split(fa_path: &Path, n: usize, out_dir: &Path, seed: u64) -> Result<()> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    let mut writers = Vec::with_capacity(n);
    for i in 1..=n {
        let path = out_dir.join(format!("{i}.fa"));
        let file = File::create(&path)
            .with_context(|| format!("failed to create {}", path.display()))?;
        writers.push(BufWriter::new(file));
    }

    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(
        File::open(fa_path).with_context(|| format!("failed to open {}", fa_path.display()))?,
    )?;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut order: Vec<usize> = (0..n).collect();

    for i in 1..=index.size {
        let j = i % n;
        if j == 0 {
            order.shuffle(&mut rng);
        }
        let seq = index.get(i)?;
        write!(&mut writers[order[j]], "{seq}")?;
    }

    for mut w in writers {
        w.flush()?;
    }

    Ok(())
}

/// Write a copy of `fa_path` with every sequence reversed. Reversed sequences
/// keep the composition of the original but destroy its homology, which is what
/// makes them usable as decoys for calibrating score cutoffs.
pub fn reverse(fa_path: &Path, out_path: &Path) -> Result<()> {
    if let Some(dir) = out_path.parent() {
        std::fs::create_dir_all(dir)?;
    }

    let mut out = BufWriter::new(
        File::create(out_path)
            .with_context(|| format!("failed to create {}", out_path.display()))?,
    );

    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(
        File::open(fa_path).with_context(|| format!("failed to open {}", fa_path.display()))?,
    )?;

    for i in 1..=index.size {
        let mut rec = index.get_record(i)?;
        rec.reverse();
        writeln!(out, "{rec}")?;
    }

    out.flush()?;
    Ok(())
}

/// Write the first `n` records of a fasta to `out_path`.
pub fn sample_to(fa_path: &Path, n: usize, out_path: &Path) -> Result<()> {
    if let Some(dir) = out_path.parent() {
        std::fs::create_dir_all(dir)?;
    }

    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(
        File::open(fa_path).with_context(|| format!("failed to open {}", fa_path.display()))?,
    )?;
    let mut out = BufWriter::new(File::create(out_path)?);

    for i in 1..=n.min(index.size) {
        write!(out, "{}", index.get(i)?)?;
    }

    out.flush()?;
    Ok(())
}
