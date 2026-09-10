//! Where this benchmark's pairs live, and where a run writes.
//!
//! Unlike the other two benchmarks there is nothing to build: the six pairs are
//! checked into `data/long-seqs/` and reached through symlinks, so `inputs/` is
//! the same on a fresh clone as it is here.
//!
//! ```text
//! inputs/query/<n>.query.fa    what is being searched with
//! inputs/target/<n>.target.fa  what it is searched against
//! outputs/                     manifest.tbl, results/, tmp/, figures/
//! ```
//!
//! A pair is named by its number, which travels through `manifest.tbl` as the
//! shard field -- one query against one target is what this benchmark measures,
//! so the pair is the unit on both sides.

use std::path::PathBuf;

use anyhow::{Context, bail};

/// This crate's directory, fixed at compile time.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn inputs() -> PathBuf {
    dir().join("inputs")
}

/// Where the run writes. No name under it: there is one pipeline and no size
/// axis, so there would be nothing to tell one output directory from another.
pub fn outputs() -> PathBuf {
    dir().join("outputs")
}

/// Where the scratch goes, one directory per thing that makes any.
///
/// Beside `outputs/` rather than inside it, so what a run produced and what it
/// merely needed on the way are not the same tree.
pub fn tmp() -> PathBuf {
    dir().join("tmp")
}

pub fn query(pair: &str) -> PathBuf {
    inputs().join(format!("query/{pair}.query.fa"))
}

pub fn target(pair: &str) -> PathBuf {
    inputs().join(format!("target/{pair}.target.fa"))
}

/// The pairs that are checked in, in numeric order.
///
/// Read off the disk rather than counted to, so adding a seventh pair is a
/// matter of dropping two files in `data/long-seqs/` -- and so a pair missing
/// half of itself is caught here by name instead of as a failed search.
pub fn pairs() -> anyhow::Result<Vec<String>> {
    let dir = inputs().join("query");
    let entries =
        std::fs::read_dir(&dir).with_context(|| format!("failed to read {}", dir.display()))?;

    // the stem of `1.query.fa` is `1.query`, so the number comes off the front
    // rather than out of the stem
    let mut out: Vec<usize> = entries
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .to_str()?
                .strip_suffix(".query.fa")?
                .parse()
                .ok()
        })
        .collect();

    out.sort_unstable();

    if out.is_empty() {
        bail!("no <n>.query.fa files in {}", dir.display());
    }

    let pairs: Vec<String> = out.into_iter().map(|n| n.to_string()).collect();

    for pair in &pairs {
        let target = target(pair);
        if !target.is_file() {
            bail!("query {pair} has no target at {}", target.display());
        }
    }

    Ok(pairs)
}
