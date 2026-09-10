//! Removing what a benchmark generated, once someone has looked at what that
//! is.
//!
//! Every benchmark keeps its generated directories at its crate root, and
//! removing them is the same job three times over: measure, show, ask, delete.
//! What differs is only which directories those are, so a benchmark hands
//! over a list and this does the rest.
//!
//! Nothing here follows a symlink. long-seqs' `inputs/` are links into
//! `data/`, and a walk that followed them would measure, and then invite
//! removing, something else entirely.

use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;

/// One directory that is about to go, and what it holds.
pub struct Target {
    /// What to call it, which is its name under the crate root.
    pub name: &'static str,
    pub path: PathBuf,
    pub files: u64,
    pub bytes: u64,
}

/// Measure what these directories hold, show it, ask, and remove them.
///
/// The ones that are not there are left out rather than reported as empty: a
/// benchmark that was never run and one that was already cleaned look the
/// same from here, and both mean there is nothing to do.
pub fn run(root: &Path, paths: &[(&'static str, PathBuf)]) -> anyhow::Result<()> {
    let mut targets = Vec::new();
    for (name, path) in paths {
        if !path.is_dir() {
            continue;
        }

        // a root that is a link belongs to whoever pointed it somewhere, and
        // what is on the other side is theirs rather than this benchmark's
        if path.is_symlink() {
            eprintln!("warning: {} is a symlink, leaving it be", path.display());
            continue;
        }

        let (files, bytes) = measure(path)?;
        targets.push(Target {
            name,
            path: path.clone(),
            files,
            bytes,
        });
    }

    if targets.is_empty() {
        println!("nothing to remove in {}", root.display());
        return Ok(());
    }

    let files: u64 = targets.iter().map(|t| t.files).sum();
    let bytes: u64 = targets.iter().map(|t| t.bytes).sum();

    println!("in {}", root.display());
    for target in &targets {
        println!(
            "  {:<10} {:>7} {:<6} {:>9}",
            target.name,
            target.files,
            plural(target.files),
            size(target.bytes)
        );
    }

    if !confirm(files, bytes)? {
        println!("left alone");
        return Ok(());
    }

    for target in &targets {
        std::fs::remove_dir_all(&target.path)
            .with_context(|| format!("failed to remove {}", target.path.display()))?;
    }

    println!("removed {} {}, {}", files, plural(files), size(bytes));
    Ok(())
}

/// How many files a directory holds, and how many bytes they take.
//
// `DirEntry::metadata` reports the link rather than its target, so a symlink
// counts as the one small file it is and is never descended into
fn measure(path: &Path) -> anyhow::Result<(u64, u64)> {
    let (mut files, mut bytes) = (0, 0);
    let mut todo = vec![path.to_path_buf()];

    while let Some(dir) = todo.pop() {
        let entries =
            std::fs::read_dir(&dir).with_context(|| format!("failed to read {}", dir.display()))?;

        for entry in entries {
            let entry = entry?;
            let meta = entry.metadata()?;

            match meta.is_dir() {
                true => todo.push(entry.path()),
                false => {
                    files += 1;
                    bytes += meta.len();
                }
            }
        }
    }

    Ok((files, bytes))
}

/// Whether to go ahead. Anything but `y` is a no, and so is a closed stdin:
/// nothing should be deleted because a script piped this somewhere.
fn confirm(files: u64, bytes: u64) -> anyhow::Result<bool> {
    print!(
        "remove {} {}, {}? [y/N] ",
        files,
        plural(files),
        size(bytes)
    );
    std::io::stdout().flush()?;

    let mut answer = String::new();
    std::io::stdin().read_line(&mut answer)?;

    Ok(matches!(answer.trim(), "y" | "Y"))
}

/// Bytes in the units `michi` writes a manifest in, so a size reads the same
/// wherever it turns up: three significant figures and no space.
fn size(bytes: u64) -> String {
    const K: f64 = 1024.0;

    let (n, unit) = match bytes as f64 {
        n if n >= K * K * K => (n / (K * K * K), "GiB"),
        n if n >= K * K => (n / (K * K), "MiB"),
        n if n >= K => (n / K, "KiB"),
        n => return format!("{n:.0}B"),
    };

    match n {
        n if n < 10.0 => format!("{n:.2}{unit}"),
        n if n < 100.0 => format!("{n:.1}{unit}"),
        n => format!("{n:.0}{unit}"),
    }
}

fn plural(files: u64) -> &'static str {
    match files {
        1 => "file",
        _ => "files",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_carry_three_figures() {
        assert_eq!(size(0), "0B");
        assert_eq!(size(512), "512B");
        assert_eq!(size(1024), "1.00KiB");
        assert_eq!(size(5 * 1024 * 1024), "5.00MiB");
        assert_eq!(size(99 * 1024 * 1024), "99.0MiB");
        assert_eq!(size(235 * 1024 * 1024), "235MiB");
        assert_eq!(size(2 * 1024 * 1024 * 1024), "2.00GiB");
    }

    #[test]
    fn a_walk_counts_files_and_not_the_directories_holding_them() {
        let dir = std::env::temp_dir().join(format!("util-clean-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(dir.join("a/b")).unwrap();
        std::fs::write(dir.join("one"), "0123456789").unwrap();
        std::fs::write(dir.join("a/two"), "01234").unwrap();
        std::fs::write(dir.join("a/b/three"), "0").unwrap();

        assert_eq!(measure(&dir).unwrap(), (3, 16));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_symlink_is_counted_but_not_followed() {
        let dir = std::env::temp_dir().join(format!("util-clean-link-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(dir.join("real")).unwrap();
        std::fs::write(dir.join("real/big"), vec![0u8; 4096]).unwrap();
        std::fs::create_dir_all(dir.join("set")).unwrap();
        std::os::unix::fs::symlink(dir.join("real"), dir.join("set/link")).unwrap();

        // the link, not the 4096 bytes on the other side of it
        let (files, bytes) = measure(&dir.join("set")).unwrap();
        assert_eq!(files, 1);
        assert!(bytes < 4096, "followed the link: {bytes} bytes");

        std::fs::remove_dir_all(&dir).ok();
    }
}
