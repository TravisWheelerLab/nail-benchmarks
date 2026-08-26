//! Where an input set lives, and what shape it has.
//!
//! Both kinds are cut from the same two sources -- Pfam and MGnify -- and
//! differ only in how. The `fixed` kind is one query set against target shards
//! of equal size, which is what a question about recall wants: the shards are
//! units of work rather than a variable. The `ladder` kind is nested rungs on
//! both axes, which is what a question about scaling wants: each rung is a
//! measurement, and every rung is a prefix of the next one up.
//!
//! A kind names a shape; a pipeline names a question. Which kind a pipeline
//! reads is fixed when it is written, so nothing here takes a kind at runtime
//! -- a pipeline reaches into its own kind's module and gets paths back.

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};

/// A way of cutting the two sources into an input set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Fixed,
    Ladder,
}

impl Kind {
    fn dir(self) -> &'static str {
        match self {
            Kind::Fixed => "fixed",
            Kind::Ladder => "ladder",
        }
    }

    pub fn queries(self) -> PathBuf {
        crate::dir().join("queries").join(self.dir())
    }

    pub fn targets(self) -> PathBuf {
        crate::dir().join("targets").join(self.dir())
    }
}

/// One query set, and targets cut into shards of equal size.
///
/// The shards are named `1..=N`, and a run over several of them is one
/// measurement over their union -- what changed between the commands is the
/// target, not the parameterization.
pub mod fixed {
    use super::*;

    const KIND: Kind = Kind::Fixed;

    pub fn queries() -> PathBuf {
        KIND.queries()
    }

    pub fn targets() -> PathBuf {
        KIND.targets()
    }

    pub fn query_hmm() -> PathBuf {
        queries().join("query.hmm")
    }

    pub fn query_sto() -> PathBuf {
        queries().join("query.sto")
    }

    /// mmseqs names a database by the file it was written to, so the path is
    /// the stem inside the directory rather than the directory.
    pub fn query_db() -> PathBuf {
        queries().join("queryDB/queryDB")
    }

    /// The target file for one shard. Takes the shard as written rather than
    /// as a number, since it travels through the manifest as a string.
    pub fn shard(shard: &str) -> PathBuf {
        targets().join(format!("{shard}.fa"))
    }

    /// Every shard, in numeric order.
    pub fn shards() -> anyhow::Result<Vec<(usize, PathBuf)>> {
        super::shards(&targets())
    }
}

/// Nested ladders on both axes: the query set at one rung is a prefix of the
/// next one up, and so is the target set.
///
/// That is what makes the grid a surface rather than a pile of unrelated
/// points. A rung is named by what was asked for -- 1000 sequences, 100
/// families -- and what it actually came to in residues is in `sizes.tbl`.
pub mod ladder {
    use super::*;

    const KIND: Kind = Kind::Ladder;

    pub fn queries() -> PathBuf {
        KIND.queries()
    }

    pub fn targets() -> PathBuf {
        KIND.targets()
    }

    /// One query rung's directory, which holds everything the tools need to
    /// search that many families.
    pub fn query(rung: usize) -> PathBuf {
        queries().join(rung.to_string())
    }

    pub fn query_hmm(rung: usize) -> PathBuf {
        query(rung).join("query.hmm")
    }

    pub fn query_sto(rung: usize) -> PathBuf {
        query(rung).join("query.sto")
    }

    pub fn query_db(rung: usize) -> PathBuf {
        query(rung).join("queryDB/queryDB")
    }

    pub fn target(rung: usize) -> PathBuf {
        targets().join(format!("{rung}.fa"))
    }

    /// The query rungs a build left behind, read off the directory rather than
    /// taken from arguments, so a run always matches the set it was pointed at.
    pub fn query_rungs() -> anyhow::Result<Vec<usize>> {
        Ok(numbered(&queries(), None)?
            .into_iter()
            .map(|(rung, _)| rung)
            .collect())
    }

    pub fn target_rungs() -> anyhow::Result<Vec<usize>> {
        Ok(numbered(&targets(), Some("fa"))?
            .into_iter()
            .map(|(rung, _)| rung)
            .collect())
    }

    /// What each rung actually came to, since sequences and models are not
    /// uniform amounts of work and residues is the honest axis to plot
    /// against. Written by `build ladder`, one file per axis.
    pub fn sizes(dir: &Path) -> PathBuf {
        dir.join("sizes.tbl")
    }
}

/// Shard files named `<n>.fa` in any directory, in numeric order.
///
/// Not every set of shards is an input set: the calibration reverses the
/// targets and shards the decoys it recruits, and reads those back the same
/// way.
pub fn shards(dir: &Path) -> anyhow::Result<Vec<(usize, PathBuf)>> {
    numbered(dir, Some("fa"))
}

/// Entries named by a number, in numeric order.
///
/// With an extension they are files named `<n>.<ext>`; without one they are
/// directories named `<n>`. The number has to come off the stem as a whole:
/// run names elsewhere embed floating point parameters, so splitting on dots is
/// not safe in general and is not worth doing differently here.
fn numbered(dir: &Path, ext: Option<&str>) -> anyhow::Result<Vec<(usize, PathBuf)>> {
    let entries =
        std::fs::read_dir(dir).with_context(|| format!("failed to read {}", dir.display()))?;

    let mut out: Vec<(usize, PathBuf)> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| match ext {
            Some(ext) => p.extension().is_some_and(|x| x == ext),
            None => p.is_dir(),
        })
        .filter_map(|p| {
            let n = p.file_stem()?.to_str()?.parse::<usize>().ok()?;
            Some((n, p))
        })
        .collect();

    out.sort_by_key(|(n, _)| *n);

    if out.is_empty() {
        let what = match ext {
            Some(ext) => format!("files named <n>.{ext}"),
            None => "directories named <n>".to_string(),
        };
        bail!("no {what} in {}; run `mgy build` first", dir.display());
    }

    Ok(out)
}
