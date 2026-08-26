//! Where an input set lives, and what it holds.
//!
//! One directory per set, named by size, holding both axes and the truth table
//! together: `build` assembles queries and targets in one pass from the same
//! pairs, and `benchmark.tbl` -- which pair is which, and at what identity --
//! belongs to neither side.
//!
//! The profmark split sits outside that, at the crate root. It is expensive,
//! it depends only on Pfam and the split parameters, and every size is drawn
//! from the same one.

use std::path::PathBuf;

use anyhow::{Context, bail};

/// This crate's directory, fixed at compile time.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The train/test split every size is drawn from.
pub fn profmark() -> PathBuf {
    dir().join("profmark")
}

pub fn profmark_query() -> PathBuf {
    profmark().join("query.sto")
}

pub fn profmark_target() -> PathBuf {
    profmark().join("target.sto")
}

/// Where the pipelines write. One directory per input set, since there is one
/// pipeline and several sizes.
pub fn runs() -> PathBuf {
    dir().join("runs")
}

/// One assembled benchmark: the queries, the targets they are hidden in, and
/// the record of which pair is which.
pub struct Inputs {
    pub size: String,
}

impl Inputs {
    pub fn new(size: &str) -> Inputs {
        Inputs {
            size: size.to_string(),
        }
    }

    pub fn dir(&self) -> PathBuf {
        dir().join("inputs").join(&self.size)
    }

    /// The profiles, built by hmmbuild from [`Inputs::query_sto`].
    pub fn query_hmm(&self) -> PathBuf {
        self.dir().join("query.hmm")
    }

    /// The query sequences, one per pair, for the tools that take sequences.
    pub fn query_fa(&self) -> PathBuf {
        self.dir().join("query.fa")
    }

    /// The alignments the profiles were built from, which mmseqs also needs.
    pub fn query_sto(&self) -> PathBuf {
        self.dir().join("query.sto")
    }

    /// One consensus sequence per profile, for asking what a profile is worth
    /// against a tool that cannot read one.
    pub fn query_cons(&self) -> PathBuf {
        self.dir().join("query.cons.fa")
    }

    /// One aligned fasta per family. psiblast takes an alignment at a time and
    /// will not read stockholm.
    pub fn afa(&self) -> PathBuf {
        self.dir().join("afa")
    }

    /// The true targets and the decoys they are hidden among.
    pub fn target_fa(&self) -> PathBuf {
        self.dir().join("target.fa")
    }

    /// Which pair is which, and at what identity. This is the benchmark's
    /// notion of truth -- there is no calibration here, and no tool is the
    /// reference.
    pub fn benchmark_tbl(&self) -> PathBuf {
        self.dir().join("benchmark.tbl")
    }

    /// Where a run over this set writes.
    pub fn run_dir(&self) -> PathBuf {
        runs().join(&self.size)
    }

    pub fn exists(&self) -> bool {
        self.dir().is_dir()
    }

    /// Every family's alignment, in name order.
    ///
    /// Sorted so a psiblast sweep runs the families in the same order every
    /// time, which is what makes two runs' wall times comparable.
    pub fn afa_files(&self) -> anyhow::Result<Vec<PathBuf>> {
        let dir = self.afa();
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
}
