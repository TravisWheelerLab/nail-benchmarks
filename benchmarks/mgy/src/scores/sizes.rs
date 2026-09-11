//! How big each target shard is, written down beside the shards themselves.
//!
//! Counting a shard means reading it. At a thousand shards of three gigabytes
//! that is the whole benchmark read a second time to fill in a metadata line,
//! so `build` writes down what it dealt and `parse` reads that.
//!
//! ```text
//! # shard count residues bytes
//! # ----- ----- -------- -----
//!   1     100000 34129933 35236472
//! ```
//!
//! A set built before this file existed has none, so a missing one is a
//! warning and a pass over the shards rather than an error.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use rayon::prelude::*;

use util::tbl;

use super::Size;

/// What the file is called, in the directory holding the shards.
pub const FILE: &str = "sizes.tbl";

pub fn path(targets: &Path) -> PathBuf {
    targets.join(FILE)
}

/// The shard a `<shard>.fa` name stands for, or the one target a pipeline with
/// no shard axis searched.
pub fn shard_path(targets: &Path, shard: &str) -> PathBuf {
    match shard.is_empty() {
        true => targets.join("target.fa"),
        false => targets.join(format!("{shard}.fa")),
    }
}

/// The size of every shard a run covered, in the order it covered them.
pub fn of(targets: &Path, shards: &[String]) -> anyhow::Result<Vec<(String, Size)>> {
    let file = path(targets);
    let mut known = match file.is_file() {
        true => read(&file)?,
        false => BTreeMap::new(),
    };

    let missing: Vec<&String> = shards
        .iter()
        .filter(|shard| !known.contains_key(*shard))
        .collect();

    if !missing.is_empty() {
        eprintln!(
            "warning: {} covers {} of the {} shards searched; reading the rest",
            file.display(),
            shards.len() - missing.len(),
            shards.len()
        );

        let shards: Vec<String> = missing.into_iter().cloned().collect();
        known.extend(measure(targets, &shards)?);
    }

    shards
        .iter()
        .map(|shard| {
            let size = known
                .get(shard)
                .with_context(|| format!("no size for shard {shard:?}"))?;

            Ok((shard.clone(), *size))
        })
        .collect()
}

/// Read every shard and count what is in it, several at a time.
pub fn measure(targets: &Path, shards: &[String]) -> anyhow::Result<Vec<(String, Size)>> {
    shards
        .par_iter()
        .map(|shard| {
            let size = fasta(&shard_path(targets, shard))?;
            Ok((shard.clone(), size))
        })
        .collect()
}

/// What every shard in a directory came to, keyed by shard.
pub fn read(path: &Path) -> anyhow::Result<BTreeMap<String, Size>> {
    let rows = tbl::read(path)?;
    let mut out = BTreeMap::new();

    for cells in &rows.cells {
        let cell = |key: &str| -> anyhow::Result<&str> {
            cells
                .get(key)
                .map(String::as_str)
                .with_context(|| format!("{} has no {key} column", path.display()))
        };

        let number = |key: &str| -> anyhow::Result<u64> {
            let text = cell(key)?;
            text.parse()
                .with_context(|| format!("bad {key} {text:?} in {}", path.display()))
        };

        let shard = match cell("shard")? {
            "-" => String::new(),
            shard => shard.to_string(),
        };

        out.insert(
            shard,
            Size {
                count: number("count")? as usize,
                residues: number("residues")?,
                bytes: number("bytes")?,
            },
        );
    }

    if out.is_empty() {
        bail!("{} has no rows", path.display());
    }

    Ok(out)
}

/// Write what a deal came to, beside the shards it wrote.
pub fn write(targets: &Path, rows: &[(String, Size)]) -> anyhow::Result<()> {
    let headers = ["shard", "count", "residues", "bytes"]
        .map(str::to_string)
        .to_vec();

    let cells: Vec<Vec<String>> = rows
        .iter()
        .map(|(shard, size)| {
            let shard = match shard.is_empty() {
                true => "-",
                false => shard,
            };

            vec![
                shard.to_string(),
                size.count.to_string(),
                size.residues.to_string(),
                size.bytes.to_string(),
            ]
        })
        .collect();

    tbl::write(
        &path(targets),
        tbl::Table {
            meta: "",
            headers: &headers,
            rows: &cells,
            ragged_last: false,
        },
    )
}

/// Records and residues in a fasta: everything that isn't a header line or
/// whitespace.
///
/// Residues is the honest unit, since neither families nor sequences are
/// uniform amounts of work.
pub fn fasta(path: &Path) -> anyhow::Result<Size> {
    let bytes = std::fs::metadata(path)
        .with_context(|| format!("failed to stat {}", path.display()))?
        .len();

    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut buf = [0u8; 1 << 16];
    let mut count = 0usize;
    let mut residues = 0u64;
    let mut in_header = false;

    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }

        for &b in &buf[..n] {
            match b {
                b'>' => {
                    in_header = true;
                    count += 1;
                }
                b'\n' => in_header = false,
                _ if in_header => {}
                _ if b.is_ascii_whitespace() => {}
                _ => residues += 1,
            }
        }
    }

    Ok(Size {
        count,
        residues,
        bytes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("mgy-sizes-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_fasta(dir: &Path, name: &str, text: &str) {
        std::fs::write(dir.join(name), text).unwrap();
    }

    #[test]
    fn sizes_survive_a_round_trip() {
        let dir = tmp("round-trip");

        let rows = vec![
            (
                "1".to_string(),
                Size {
                    count: 3,
                    residues: 120,
                    bytes: 148,
                },
            ),
            (
                String::new(),
                Size {
                    count: 1,
                    residues: 7,
                    bytes: 20,
                },
            ),
        ];

        write(&dir, &rows).unwrap();
        let back = read(&path(&dir)).unwrap();

        assert_eq!(back.len(), 2);
        assert_eq!(back["1"].residues, 120);
        assert_eq!(back[""].count, 1);
    }

    #[test]
    fn a_shard_with_no_row_is_read_off_the_file() {
        let dir = tmp("fallback");

        write_fasta(&dir, "1.fa", ">a\nACDEF\nGH\n>b\nIKLM\n");
        write_fasta(&dir, "2.fa", ">c\nAAA\n");

        let measured = of(&dir, &["1".to_string(), "2".to_string()]).unwrap();
        assert_eq!(measured[0].1.count, 2);
        assert_eq!(measured[0].1.residues, 11);
        assert_eq!(measured[1].1.residues, 3);

        // written down, the second pass reads rather than measures
        write(&dir, &measured).unwrap();
        let known = of(&dir, &["1".to_string(), "2".to_string()]).unwrap();

        assert_eq!(known[0].1.count, measured[0].1.count);
        assert_eq!(known[0].1.residues, measured[0].1.residues);
        assert_eq!(known[1].1.bytes, measured[1].1.bytes);
    }
}
