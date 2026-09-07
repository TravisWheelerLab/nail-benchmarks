//! Cutting a query set into files a batch of jobs can search in parallel.

use std::borrow::Borrow;
use std::cmp::Reverse;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use libsail::collection::{Indexable, Iterable};
use libsail::seq::fasta::IndexedFasta;
use libsail::seq::p7hmm::IndexedHmm;

/// Which format is being cut up, chosen by the caller at runtime.
#[derive(Clone, Copy, Debug)]
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

/// Write `n` splits of `path` into `out_dir`, returning their paths.
///
/// Empty parts are skipped, so fewer files than `n` come back when there are
/// fewer records than requested splits.
pub fn write_splits(
    path: impl AsRef<Path>,
    kind: Kind,
    n: usize,
    out_dir: impl AsRef<Path>,
) -> anyhow::Result<Vec<PathBuf>> {
    let path = path.as_ref();
    let out_dir = out_dir.as_ref();

    if n == 0 {
        bail!("cannot split {} into 0 parts", path.display());
    }

    let opened = || format!("failed to index {}", path.display());
    let empty = || format!("no records found in {}", path.display());

    match kind {
        Kind::Hmm => {
            let c = IndexedHmm::open(path).with_context(opened)?;
            ensure!(!c.is_empty(), empty());
            deal(c, n, out_dir, kind.extension(), |rec| rec.header.leng)
        }
        Kind::Fasta => {
            let c = IndexedFasta::open(path).with_context(opened)?;
            ensure!(!c.is_empty(), empty());
            deal(c, n, out_dir, kind.extension(), |rec| rec.seq.len())
        }
    }
}

/// Deal `c` into `n` files, `<i>.<ext>` for `i` in `0..n`, heaviest record
/// first onto whichever part is lightest so far.
///
/// Weight rather than count: a model's length and a sequence's residues are
/// what a search against it costs, and the parts go to jobs whose slowest one
/// is the wall time this repo reports. Round robin over the same order leaves
/// 0.5% of the work misplaced at 8 parts and 14% at 256; this leaves 0.03%.
fn deal<C>(
    c: C,
    n: usize,
    out_dir: &Path,
    ext: &str,
    weight: impl Fn(&C::Record) -> usize,
) -> anyhow::Result<Vec<PathBuf>>
where
    C: Indexable,
    C::Record: Display,
{
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    // descending, so the walk below is longest-processing-time order and each
    // record's own weight is in hand when its part is chosen
    let sorted = c.sort_by_key(|rec| Reverse(weight(rec)));
    let parts = n.min(sorted.len());

    let mut writers = Vec::with_capacity(parts);
    let mut written = Vec::with_capacity(parts);

    for i in 0..parts {
        let path = out_dir.join(format!("{i}.{ext}"));
        let file = File::create(&path)
            .with_context(|| format!("failed to create {}", path.display()))?;

        writers.push(BufWriter::new(file));
        written.push(path);
    }

    let mut loads = vec![0usize; parts];

    for rec in sorted.iter() {
        let rec = rec.borrow();

        // ties to the earlier part, so a given input deals the same way twice
        let lightest = loads
            .iter()
            .enumerate()
            .min_by_key(|(i, load)| (**load, *i))
            .map(|(i, _)| i)
            .expect("at least one part");

        loads[lightest] += weight(rec);
        write!(writers[lightest], "{rec}")?;
    }

    for mut w in writers {
        w.flush()?;
    }

    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profile;
    use libsail::collection::Iterable;

    // one directory per test: they run in parallel, and a shared one gets
    // removed out from under its neighbours
    fn tmp(name: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "bench-split-{}-{}",
            std::process::id(),
            name.replace('.', "_")
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, body).unwrap();
        path
    }

    fn names_in(path: &Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter_map(|l| l.strip_prefix("NAME").map(|r| r.trim().to_string()))
            .collect()
    }

    #[test]
    fn weighs_a_profile_by_its_model_length() {
        let path = tmp("a.hmm", &format!("{}{}", profile("one", 3), profile("two", 7)));
        let c = IndexedHmm::open(&path).unwrap();

        let lengs: Vec<usize> = c.iter().map(|r| r.header.leng).collect();
        assert_eq!(lengs, vec![3, 7]);

        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }

    #[test]
    fn weighs_a_sequence_by_its_residues() {
        let path = tmp("a.fa", ">one\nAAAA\nCC\n>two\nDDD\n");
        let c = IndexedFasta::open(&path).unwrap();

        let lens: Vec<usize> = c.iter().map(|r| r.seq.len()).collect();
        assert_eq!(lens, vec![6, 3]);

        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }

    #[test]
    fn splits_round_trip_every_record() {
        let body: String = ["a", "b", "c"]
            .iter()
            .enumerate()
            .map(|(i, n)| profile(n, (i + 1) * 10))
            .collect();

        let path = tmp("b.hmm", &body);
        let out = path.parent().unwrap().join("splits");
        let paths = write_splits(&path, Kind::Hmm, 2, &out).unwrap();

        assert_eq!(paths.len(), 2);

        let mut names: Vec<String> = paths.iter().flat_map(|p| names_in(p)).collect();
        names.sort();
        assert_eq!(names, vec!["a", "b", "c"]);

        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }

    #[test]
    fn a_split_is_a_profile_a_search_can_read_back() {
        let path = tmp("d.hmm", &profile("solo", 4));
        let out = path.parent().unwrap().join("splits-d");
        let paths = write_splits(&path, Kind::Hmm, 1, &out).unwrap();

        let back = IndexedHmm::open(&paths[0]).unwrap();
        assert_eq!(back.len(), 1);

        let rec = back.cloned(0).unwrap();
        assert_eq!(rec.header.name, "solo");
        assert_eq!(rec.header.leng, 4);
        assert_eq!(rec.model.match_emissions.len(), 4 * 2);

        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }

    #[test]
    fn dealing_spreads_weight_evenly() {
        // 6 profiles of 10..=60 nodes: 210 nodes over 3 parts balances at 70
        let body: String = (1..=6).map(|i| profile(&format!("p{i}"), i * 10)).collect();

        let path = tmp("e.hmm", &body);
        let out = path.parent().unwrap().join("splits-e");
        let paths = write_splits(&path, Kind::Hmm, 3, &out).unwrap();

        let loads: Vec<usize> = paths
            .iter()
            .map(|p| {
                IndexedHmm::open(p)
                    .unwrap()
                    .iter()
                    .map(|r| r.header.leng)
                    .sum()
            })
            .collect();

        assert_eq!(loads, vec![70, 70, 70]);

        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }

    #[test]
    fn fewer_records_than_splits_yields_fewer_files() {
        let path = tmp("c.fa", ">only\nAAAA\n");
        let out = path.parent().unwrap().join("splits-c");
        let paths = write_splits(&path, Kind::Fasta, 4, &out).unwrap();

        assert_eq!(paths.len(), 1);
        assert_eq!(names_in(&paths[0]).len(), 0, "a fasta has no NAME lines");

        std::fs::remove_dir_all(path.parent().unwrap()).ok();
    }
}
