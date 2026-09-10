//! Where the benchmark lives, and what it holds.
//!
//! One directory, holding both axes and the truth table together: `build`
//! assembles queries and targets in one pass from the same pairs, and
//! `benchmark.tbl` -- which pair is which, and at what identity -- belongs to
//! neither side.
//!
//! The profmark split sits outside it, at the crate root. It is expensive and
//! depends only on Pfam and the split parameters, so rebuilding the benchmark
//! draws from the same split rather than making a new one.

use std::path::PathBuf;

use anyhow::{Context, bail};

/// This crate's directory, fixed at compile time.
pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The train/test split the benchmark is drawn from.
pub fn profmark() -> PathBuf {
    root().join("profmark")
}

pub fn profmark_query() -> PathBuf {
    profmark().join("query.sto")
}

pub fn profmark_target() -> PathBuf {
    profmark().join("target.sto")
}

/// The assembled benchmark: the queries, the targets they are hidden in, and
/// the record of which pair is which.
pub fn dir() -> PathBuf {
    root().join("inputs")
}

/// Where a run writes.
pub fn outputs() -> PathBuf {
    root().join("outputs")
}

/// Where the scratch goes, one directory per thing that makes any.
///
/// Beside `outputs/` rather than inside it, so what a run produced and what it
/// merely needed on the way are not the same tree.
pub fn tmp() -> PathBuf {
    root().join("tmp")
}

pub fn exists() -> bool {
    dir().is_dir()
}

/// The profiles, built by hmmbuild from [`query_sto`].
pub fn query_hmm() -> PathBuf {
    dir().join("query.hmm")
}

/// The query sequences, one per pair, for the tools that take sequences.
pub fn query_fa() -> PathBuf {
    dir().join("query.fa")
}

/// The alignments the profiles were built from, which mmseqs also needs.
pub fn query_sto() -> PathBuf {
    dir().join("query.sto")
}

/// One consensus sequence per profile, written by `build`.
//
// nothing under `run` searches against it: the modes are prf and seq, so this
// is written and never read
pub fn query_cons() -> PathBuf {
    dir().join("query.cons.fa")
}

/// One aligned fasta per family. psiblast takes an alignment at a time and
/// will not read stockholm.
pub fn afa() -> PathBuf {
    dir().join("afa")
}

/// The true targets and the decoys they are hidden among.
pub fn target_fa() -> PathBuf {
    dir().join("target.fa")
}

/// Which pair is which, and at what identity. This is the benchmark's notion
/// of truth -- there is no calibration here, and no tool is the reference.
pub fn benchmark_tbl() -> PathBuf {
    dir().join("benchmark.tbl")
}

/// Every family's alignment, in name order.
///
/// Sorted so a psiblast sweep runs the families in the same order every time,
/// which is what makes two runs' wall times comparable.
pub fn afa_files() -> anyhow::Result<Vec<PathBuf>> {
    let dir = afa();
    let entries =
        std::fs::read_dir(&dir).with_context(|| format!("failed to read {}", dir.display()))?;

    let mut out: Vec<PathBuf> = entries
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "afa"))
        .collect();

    out.sort();

    if out.is_empty() {
        bail!("no .afa files in {}", dir.display());
    }

    Ok(out)
}
