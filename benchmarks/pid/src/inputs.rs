//! The one unit of a profmark set, with its artifacts resolved.
//!
//! `build-set --in pid-toy|pid-real` writes the set; this reads it. The whole
//! benchmark is a single unit — every query against one target file — so there
//! is one row, and what tells a true pair from a decoy lives in the file the
//! `truth` column names rather than in a column of its own.

use std::path::{Path, PathBuf};

use anyhow::Context;

use util::set::{Set, shape};

pub struct Inputs {
    /// The profiles, built by hmmbuild from [`Inputs::query_sto`].
    pub query_hmm: PathBuf,
    /// The query sequences, one per pair, for the tools that take sequences.
    pub query_fa: PathBuf,
    /// The alignments the profiles were built from, which mmseqs also needs.
    pub query_sto: PathBuf,
    /// One consensus sequence per profile.
    //
    // nothing under `run` searches against it: the modes are prf and seq, so
    // the recipe writes it and nothing here reads it
    #[allow(dead_code)]
    pub query_cons: PathBuf,
    /// One aligned fasta per family. psiblast takes an alignment at a time and
    /// will not read stockholm.
    pub afa: PathBuf,
    /// The true targets and the decoys they are hidden among.
    pub target_fa: PathBuf,
    /// Which pair is which, and at what identity. The benchmark's own notion
    /// of truth: there is no calibration here and no tool is the reference.
    pub truth: PathBuf,
}

impl Inputs {
    pub fn open(set_dir: &Path) -> anyhow::Result<Inputs> {
        let set = Set::load_as(set_dir, &shape::PROFMARK)?;

        let unit = set
            .units()
            .next()
            .with_context(|| format!("{} names no units", set_dir.display()))?;

        let query_hmm = unit.query_hmm()?;
        // beside the profiles rather than in the manifest: both are written by
        // the recipe and neither is a representation a search asks for
        let queries = query_hmm
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf();

        Ok(Inputs {
            query_fa: unit.query_fa()?,
            query_sto: unit.query_sto()?,
            query_cons: queries.join("query.cons.fa"),
            afa: queries.join("afa"),
            target_fa: unit.target()?,
            truth: set_dir.join(unit.need("truth")?),
            query_hmm,
        })
    }

    /// Every family's alignment, in name order.
    ///
    /// Sorted so a psiblast sweep runs the families in the same order every
    /// time, which is what makes two runs' wall times comparable.
    pub fn afa_files(&self) -> anyhow::Result<Vec<PathBuf>> {
        let mut out: Vec<PathBuf> = std::fs::read_dir(&self.afa)
            .with_context(|| format!("failed to read {}", self.afa.display()))?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "afa"))
            .collect();

        out.sort();
        Ok(out)
    }
}
