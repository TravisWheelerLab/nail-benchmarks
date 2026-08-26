use std::path::{Path, PathBuf};

use anyhow::{bail, Context};

pub use crate::dir;

pub fn shards(dir: &Path) -> anyhow::Result<Vec<(usize, PathBuf)>> {
    let entries =
        std::fs::read_dir(dir).with_context(|| format!("failed to read {}", dir.display()))?;

    let mut out: Vec<(usize, PathBuf)> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .filter_map(|p| {
            let i = p.file_stem()?.to_str()?.parse::<usize>().ok()?;
            Some((i, p))
        })
        .collect();

    out.sort_by_key(|(i, _)| *i);

    if out.is_empty() {
        bail!("no shards named <n>.fa in {}", dir.display());
    }

    Ok(out)
}
