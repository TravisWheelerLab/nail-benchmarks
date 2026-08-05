use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

/// One indexed record (an HMM block or a FASTA entry) as a byte range plus the
/// weight used to balance splits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Record {
    pub start: u64,
    pub end: u64,
    pub weight: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Hmm,
    Fasta,
}

impl Kind {
    pub fn extension(&self) -> &'static str {
        match self {
            Kind::Hmm => "hmm",
            Kind::Fasta => "fa",
        }
    }
}

/// Index records without retaining file contents. Query sets get large (the
/// mgnify benchmark's query.hmm is the full 1.6GB Pfam), so this pass keeps
/// only byte offsets and writes splits by seeking back over the original.
pub fn index(path: impl AsRef<Path>, kind: Kind) -> Result<Vec<Record>> {
    let path = path.as_ref();
    let file =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut line = Vec::new();
    let mut pos: u64 = 0;

    let mut start: u64 = 0;
    let mut weight: u64 = 0;
    let mut open = false;

    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("failed reading {}", path.display()))?;
        if n == 0 {
            break;
        }

        let end = pos + n as u64;

        match kind {
            Kind::Hmm => {
                if line.starts_with(b"//") {
                    records.push(Record {
                        start,
                        end,
                        weight,
                    });
                    start = end;
                    weight = 0;
                    open = false;
                } else {
                    open = true;
                    if line.starts_with(b"LENG") {
                        weight = parse_field(&line, 1).unwrap_or(0);
                    }
                }
            }
            Kind::Fasta => {
                if line.starts_with(b">") {
                    if open {
                        records.push(Record {
                            start,
                            end: pos,
                            weight,
                        });
                    }
                    start = pos;
                    weight = 0;
                    open = true;
                } else if open {
                    weight += line
                        .iter()
                        .filter(|b| !b.is_ascii_whitespace())
                        .count() as u64;
                }
            }
        }

        pos = end;
    }

    // trailing fasta record has no following '>' to close it
    if kind == Kind::Fasta && open {
        records.push(Record {
            start,
            end: pos,
            weight,
        });
    }

    if records.is_empty() {
        bail!("no records found in {}", path.display());
    }

    Ok(records)
}

fn parse_field(line: &[u8], idx: usize) -> Option<u64> {
    std::str::from_utf8(line)
        .ok()?
        .split_whitespace()
        .nth(idx)?
        .parse()
        .ok()
}

/// Distribute records into `n` bins, heaviest first onto the lightest bin
/// (longest-processing-time). Replaces the sort-then-round-robin deal in
/// util/scripts/common.py, which balanced the same quantity slightly worse.
pub fn balance(records: &[Record], n: usize) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_by(|&a, &b| {
        records[b]
            .weight
            .cmp(&records[a].weight)
            .then(a.cmp(&b))
    });

    let mut bins = vec![Vec::new(); n];
    let mut loads = vec![0u64; n];

    for idx in order {
        let lightest = loads
            .iter()
            .enumerate()
            .min_by_key(|(i, load)| (**load, *i))
            .map(|(i, _)| i)
            .expect("at least one bin");

        bins[lightest].push(idx);
        loads[lightest] += records[idx].weight;
    }

    // keep each bin in original file order for reproducible output
    for bin in &mut bins {
        bin.sort_unstable();
    }

    bins
}

/// Write `n` balanced splits of `path` into `out_dir`, returning their paths.
/// Empty bins are skipped, so fewer files than `n` may come back when there
/// are fewer records than requested splits.
pub fn write_splits(
    path: impl AsRef<Path>,
    kind: Kind,
    n: usize,
    out_dir: impl AsRef<Path>,
) -> Result<Vec<PathBuf>> {
    let path = path.as_ref();
    let out_dir = out_dir.as_ref();

    if n == 0 {
        bail!("cannot split {} into 0 parts", path.display());
    }

    let records = index(path, kind)?;
    let bins = balance(&records, n);

    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    let mut src = File::open(path)?;
    let mut written = Vec::new();

    for (i, bin) in bins.iter().enumerate() {
        if bin.is_empty() {
            continue;
        }

        let out_path = out_dir.join(format!("{i}.{}", kind.extension()));
        let mut out = BufWriter::new(
            File::create(&out_path)
                .with_context(|| format!("failed to create {}", out_path.display()))?,
        );

        for &idx in bin {
            let rec = records[idx];
            src.seek(SeekFrom::Start(rec.start))?;
            let mut chunk = (&mut src).take(rec.end - rec.start);
            std::io::copy(&mut chunk, &mut out)?;
        }

        out.flush()?;
        written.push(out_path);
    }

    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    // one directory per test: they run in parallel, and a shared one gets
    // removed out from under its neighbours
    fn tmp(name: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "bioio-split-{}-{}",
            std::process::id(),
            name.replace('.', "_")
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        let mut f = File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        path
    }

    #[test]
    fn indexes_hmm_blocks_by_model_length() {
        let path = tmp(
            "a.hmm",
            "HMMER3/f\nNAME  one\nLENG  100\n//\nHMMER3/f\nNAME  two\nLENG  250\n//\n",
        );
        let records = index(&path, Kind::Hmm).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].weight, 100);
        assert_eq!(records[1].weight, 250);
        // blocks must tile the file with no gaps or overlap
        assert_eq!(records[0].start, 0);
        assert_eq!(records[0].end, records[1].start);
    }

    #[test]
    fn indexes_fasta_by_residue_count() {
        let path = tmp("a.fa", ">one\nAAAA\nCC\n>two\nDDD\n");
        let records = index(&path, Kind::Fasta).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].weight, 6);
        assert_eq!(records[1].weight, 3);
    }

    #[test]
    fn balance_spreads_weight_evenly() {
        let records: Vec<Record> = [10u64, 9, 8, 7, 6, 5]
            .iter()
            .map(|&w| Record {
                start: 0,
                end: 0,
                weight: w,
            })
            .collect();

        let bins = balance(&records, 3);
        assert_eq!(bins.len(), 3);

        let loads: Vec<u64> = bins
            .iter()
            .map(|b| b.iter().map(|&i| records[i].weight).sum())
            .collect();

        // 45 total across 3 bins balances exactly at 15 each
        assert_eq!(loads, vec![15, 15, 15]);
    }

    #[test]
    fn splits_round_trip_every_record() {
        let path = tmp(
            "b.hmm",
            "HMMER3/f\nNAME  a\nLENG  10\n//\nHMMER3/f\nNAME  b\nLENG  20\n//\nHMMER3/f\nNAME  c\nLENG  30\n//\n",
        );
        let out = path.parent().unwrap().join("splits");
        let paths = write_splits(&path, Kind::Hmm, 2, &out).unwrap();

        assert_eq!(paths.len(), 2);

        let mut names: Vec<String> = Vec::new();
        for p in &paths {
            for line in std::fs::read_to_string(p).unwrap().lines() {
                if let Some(rest) = line.strip_prefix("NAME") {
                    names.push(rest.trim().to_string());
                }
            }
        }
        names.sort();

        assert_eq!(names, vec!["a", "b", "c"]);
        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }

    #[test]
    fn fewer_records_than_splits_yields_fewer_files() {
        let path = tmp("c.fa", ">only\nAAAA\n");
        let out = path.parent().unwrap().join("splits-c");
        let paths = write_splits(&path, Kind::Fasta, 4, &out).unwrap();

        assert_eq!(paths.len(), 1);
        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }
}
