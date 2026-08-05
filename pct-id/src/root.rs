use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

/// Locate the repository root by walking up from a starting point until a
/// directory holding the workspace `Cargo.toml` is found.
///
/// Deliberately a runtime search rather than `env!("CARGO_MANIFEST_DIR")`:
/// benchmarks are built on one machine and run on another, and a compile-time
/// path would point somewhere that does not exist there.
pub fn find(override_root: Option<&Path>) -> Result<PathBuf> {
    if let Some(root) = override_root {
        let root = root
            .canonicalize()
            .with_context(|| format!("--root {} does not exist", root.display()))?;
        return Ok(root);
    }

    // start from the executable so an installed or copied binary still works,
    // then fall back to the working directory
    let starts = [
        std::env::current_exe().ok().and_then(|p| p.parent().map(Path::to_path_buf)),
        std::env::current_dir().ok(),
    ];

    for start in starts.into_iter().flatten() {
        if let Some(root) = walk_up(&start) {
            return Ok(root);
        }
    }

    bail!(
        "could not locate the repository root; run from inside the repo or pass --root"
    )
}

fn walk_up(start: &Path) -> Option<PathBuf> {
    for dir in start.ancestors() {
        let manifest = dir.join("Cargo.toml");
        if !manifest.exists() {
            continue;
        }
        // only the workspace manifest declares [workspace]; member manifests
        // would otherwise stop the walk one level too early
        if std::fs::read_to_string(&manifest)
            .map(|t| t.contains("[workspace]"))
            .unwrap_or(false)
        {
            return Some(dir.to_path_buf());
        }
    }
    None
}
