//! How far down the prefilter list the hits keep coming.
//!
//! An arm seeded static and unbounded aligned everything its prefilter
//! returned, so every pair nail reported has a depth: its rank in its query's
//! prefilter list, which mmseqs writes in score order. Depth is binned in
//! ranges that double, the way prog's `n_take` does, and each bin counts the
//! prefilter pairs at that depth against the seeds, the pairs the arm kept
//! over the cutoff, the pairs the ceiling kept, and the hits: kept by both.
//! The hits over the prefilter pairs is what a stopping rule is betting on.
//!
//! The prefilter list is read out of the databases nail left under
//! `results/prefilter.<arm>.<shard>/`: mmseqs' index gives each query's
//! entry, split across one data file per thread whose offsets run on as if
//! the files were one, and nail's header databases give which key is which
//! name, in the order it wrote them.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};

use crate::scores::runs::Reader as Runs;

/// Where one arm's seeding left mmseqs' databases, under a run's root.
pub fn prefilter_dir(root: &Path, arm: &str, shard: &str) -> PathBuf {
    root.join("results")
        .join(format!("prefilter.{arm}.{shard}"))
}

/// One arm's prefilter database for one shard, open for reading ranks.
struct Prefilter {
    /// Each query key's entry: where it starts in the joined data, and how
    /// long it is.
    entries: HashMap<u32, (u64, u64)>,
    /// The data files in mmseqs' order, each with where it starts in the
    /// joined data and its length.
    files: Vec<(File, u64, u64)>,
    query_key: HashMap<String, u32>,
    target_key: HashMap<String, u32>,
}

impl Prefilter {
    fn open(dir: &Path) -> anyhow::Result<Prefilter> {
        let pdb = dir.join("prefilter-db");
        let entries = index(&pdb.join("pdb.index"))?
            .into_iter()
            .map(|(key, off, len)| (key, (off, len)))
            .collect();

        // pdb.0, pdb.1, ... pdb.15: numeric order, since the offsets in the
        // index run across the files in the order mmseqs wrote them
        let mut parts: Vec<(u32, PathBuf)> = std::fs::read_dir(&pdb)
            .with_context(|| format!("failed to read {}", pdb.display()))?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter_map(|p| {
                let n: u32 = p.extension()?.to_str()?.parse().ok()?;
                (p.file_stem()? == "pdb").then_some((n, p))
            })
            .collect();
        parts.sort_unstable();
        ensure!(
            !parts.is_empty(),
            "no pdb.N data files in {}",
            pdb.display()
        );

        let mut files = Vec::with_capacity(parts.len());
        let mut start = 0u64;
        for (_, path) in parts {
            let file =
                File::open(&path).with_context(|| format!("failed to open {}", path.display()))?;
            let len = file.metadata()?.len();
            files.push((file, start, len));
            start += len;
        }

        Ok(Prefilter {
            entries,
            files,
            query_key: headers(&dir.join("query-db").join("qdb_h"))?,
            target_key: headers(&dir.join("target-db").join("tdb_h"))?,
        })
    }

    /// One query's entry, as the lines mmseqs wrote.
    fn entry(&self, key: u32) -> anyhow::Result<Vec<u8>> {
        let &(off, len) = self
            .entries
            .get(&key)
            .with_context(|| format!("no prefilter entry for query key {key}"))?;

        let mut out = vec![0u8; len as usize];
        let mut done = 0u64;
        // an entry sits inside one file, but reading across the seam
        // costs nothing and holds if that ever changes
        while done < len {
            let at = off + done;
            let (file, start, flen) = self
                .files
                .iter()
                .find(|(_, start, flen)| at >= *start && at < start + flen)
                .with_context(|| format!("offset {at} is past the prefilter data"))?;
            let take = (len - done).min(start + flen - at);
            file.read_exact_at(&mut out[done as usize..(done + take) as usize], at - start)?;
            done += take;
        }
        Ok(out)
    }

    /// Each target key in one query's list and its rank there, from 1.
    fn ranks(&self, key: u32) -> anyhow::Result<HashMap<u32, u32>> {
        let entry = self.entry(key)?;
        let mut ranks = HashMap::new();
        for (i, line) in entry.split(|&b| b == b'\n').enumerate() {
            let Some(field) = line.split(|&b| b == b'\t').next() else {
                continue;
            };
            if field.is_empty() || field == b"\0" {
                continue;
            }
            let target: u32 = std::str::from_utf8(field)?
                .parse()
                .with_context(|| format!("bad target key in prefilter entry {key}"))?;
            ranks.insert(target, i as u32 + 1);
        }
        Ok(ranks)
    }

    /// How many targets each query's list holds, by key, by reading the
    /// data once through.
    fn lengths(&self) -> anyhow::Result<Vec<(u32, u64)>> {
        let mut spans: Vec<(u64, u64, u32)> = self
            .entries
            .iter()
            .map(|(&key, &(off, len))| (off, len, key))
            .collect();
        spans.sort_unstable();

        let mut lengths = Vec::with_capacity(spans.len());
        let mut buf = vec![0u8; 8 << 20];

        for (file, start, flen) in &self.files {
            let mut reader = BufReader::with_capacity(buf.len(), file);
            let mut at = *start;
            let end = start + flen;

            while at < end {
                let n = reader.read(&mut buf)?;
                if n == 0 {
                    break;
                }
                let chunk = &buf[..n];

                // count the newlines of every span this chunk touches
                let mut i = spans.partition_point(|&(off, len, _)| off + len <= at);
                while i < spans.len() && spans[i].0 < at + n as u64 {
                    let (off, len, _) = spans[i];
                    let lo = off.max(at) - at;
                    let hi = (off + len).min(at + n as u64) - at;
                    let count = chunk[lo as usize..hi as usize]
                        .iter()
                        .filter(|&&b| b == b'\n')
                        .count() as u64;
                    if lengths.len() <= i {
                        lengths.resize(i + 1, 0);
                    }
                    lengths[i] += count;
                    i += 1;
                }

                at += n as u64;
            }
        }

        lengths.resize(spans.len(), 0);
        Ok(spans
            .iter()
            .zip(lengths)
            .map(|(&(_, _, key), len)| (key, len))
            .collect())
    }
}

/// An mmseqs index: key, offset, length per line.
fn index(path: &Path) -> anyhow::Result<Vec<(u32, u64, u64)>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let mut f = line.split_whitespace();
        let (Some(k), Some(o), Some(l)) = (f.next(), f.next(), f.next()) else {
            bail!("short line in {}: {line:?}", path.display());
        };
        out.push((k.parse()?, o.parse()?, l.parse()?));
    }
    Ok(out)
}

/// Which key nail gave each name, from a header database it wrote: one
/// header per entry, the name up to the first whitespace.
fn headers(path: &Path) -> anyhow::Result<HashMap<String, u32>> {
    let name = path.file_name().context("a header database needs a name")?;
    let idx = index(&path.with_file_name(format!("{}.index", name.to_string_lossy())))?;
    let mut data = Vec::new();
    File::open(path)
        .with_context(|| format!("failed to open {}", path.display()))?
        .read_to_end(&mut data)?;

    let mut map = HashMap::with_capacity(idx.len());
    for (key, off, len) in idx {
        let entry = &data[off as usize..(off + len) as usize];
        let name = entry
            .split(|b| b.is_ascii_whitespace() || *b == 0)
            .next()
            .unwrap_or_default();
        map.insert(String::from_utf8_lossy(name).into_owned(), key);
    }
    Ok(map)
}

// ---

/// The bin a rank falls in: the first holds ranks up to `width`, and each
/// after it doubles the ceiling.
fn bin_of(rank: u64, width: u64) -> usize {
    let mut hi = width;
    let mut k = 0;
    while rank > hi {
        hi *= 2;
        k += 1;
    }
    k
}

/// The rank range of one bin.
fn range(bin: usize, width: u64) -> (u64, u64) {
    match bin {
        0 => (1, width),
        k => ((width << (k - 1)) + 1, width << k),
    }
}

/// One (unit, run)'s counts per bin.
#[derive(Default)]
struct Tally {
    prefilter: Vec<u64>,
    seeds: Vec<u64>,
    kept: Vec<u64>,
    ceiling: Vec<u64>,
    /// Kept by this arm and by the ceiling: the pairs a stopping rule is
    /// betting on.
    hits: Vec<u64>,
    /// The ceiling's hits this arm's prefilter never returned, so they have
    /// no depth in it.
    ceiling_beyond: u64,
}

impl Tally {
    fn bump(v: &mut Vec<u64>, bin: usize) {
        if v.len() <= bin {
            v.resize(bin + 1, 0);
        }
        v[bin] += 1;
    }

    fn bins(&self) -> usize {
        [
            &self.prefilter,
            &self.seeds,
            &self.kept,
            &self.ceiling,
            &self.hits,
        ]
        .iter()
        .map(|v| v.len())
        .max()
        .unwrap_or(0)
    }
}

/// A pair one run kept, or the ceiling kept, waiting on its rank.
struct Pair {
    query: String,
    target: String,
    kept: bool,
    ceiling: bool,
}

/// Where the hits sit in the prefilter list, per unit and per run, and every
/// hit's own rank in `hits`.
pub fn depth(
    table: &Path,
    root: &Path,
    width: u64,
    out: &Path,
    hits_out: &Path,
    lists_out: &Path,
) -> anyhow::Result<()> {
    ensure!(width > 0, "--bin must be at least 1");

    let mut scores = Runs::open(table)?;
    let top = scores.meta().ceiling()?;
    let runs = scores.meta().runs.clone();

    // every pair the ceiling or a run kept, grouped by the arm and the shard
    // whose prefilter list says where it sat
    let mut pairs: HashMap<(usize, String), Vec<Pair>> = HashMap::new();
    // what the ceiling kept per unit, the denominator a sensitivity is read
    // over
    let mut hits_of: indexmap::IndexMap<String, u64> = indexmap::IndexMap::new();

    scores.each(|row| {
        let ceiling = row.row().passed(top);
        if ceiling {
            *hits_of.entry(row.row().shard().to_string()).or_default() += 1;
        }
        let query = String::from_utf8_lossy(row.row().field(0)).into_owned();
        let target = String::from_utf8_lossy(row.row().field(1)).into_owned();

        for (run, column) in runs.iter().enumerate() {
            if column.seeds.is_none() {
                continue;
            }
            let kept = row.present(run) && row.row().passed(run);
            if !kept && !ceiling {
                continue;
            }
            pairs
                .entry((run, row.row().shard().to_string()))
                .or_default()
                .push(Pair {
                    query: query.clone(),
                    target: target.clone(),
                    kept,
                    ceiling,
                });
        }
        Ok(())
    })?;

    let mut table = toil::Table::new(toil::Schema::new([
        toil::Column::new("unit"),
        toil::Column::new("run"),
        toil::Column::new("lo"),
        toil::Column::new("hi"),
        toil::Column::new("prefilter"),
        toil::Column::new("seeds"),
        toil::Column::new("kept"),
        toil::Column::new("ceiling"),
        toil::Column::new("hits"),
        toil::Column::new("frac").fixed(4),
    ]));
    table.meta("bin", [width.to_string()]);
    table.meta("ceiling", [runs[top].name.as_str()]);
    let mut hits_table = toil::Table::new(toil::Schema::new([
        toil::Column::new("unit"),
        toil::Column::new("run"),
        toil::Column::new("query"),
        toil::Column::new("target"),
        toil::Column::new("rank"),
    ]));
    let mut lists_table = toil::Table::new(toil::Schema::new([
        toil::Column::new("unit"),
        toil::Column::new("run"),
        toil::Column::new("query"),
        toil::Column::new("length"),
    ]));
    hits_table.meta("ceiling", [runs[top].name.as_str()]);
    for (unit, n) in &hits_of {
        table.meta("hits", [unit.clone(), n.to_string()]);
        hits_table.meta("hits", [unit.clone(), n.to_string()]);
    }

    let mut groups: Vec<_> = pairs.into_iter().collect();
    groups.sort_by(|a, b| (a.0.1.as_str(), a.0.0).cmp(&(b.0.1.as_str(), b.0.0)));

    for ((run, shard), mut pairs) in groups {
        let arm = runs[run].seeds.as_deref().unwrap();
        let dir = prefilter_dir(root, arm, &shard);
        let pf = Prefilter::open(&dir).with_context(|| {
            format!("failed to open the prefilter database in {}", dir.display())
        })?;

        let mut tally = Tally::default();

        // the prefilter pairs per bin come from the length of every
        // query's list: a list of L holds min(L, hi) - lo + 1 of each bin
        // the name each key was given, for the per-query list
        let name_of: HashMap<u32, &str> = pf
            .query_key
            .iter()
            .map(|(name, &key)| (key, name.as_str()))
            .collect();

        for (key, len) in pf.lengths()? {
            lists_table.row([
                toil::Cell::from(shard.as_str()),
                runs[run].name.as_str().into(),
                (*name_of.get(&key).unwrap_or(&"?")).into(),
                len.into(),
            ]);
            let mut bin = 0;
            loop {
                let (lo, hi) = range(bin, width);
                if len < lo {
                    break;
                }
                Tally::bump(&mut tally.prefilter, bin);
                tally.prefilter[bin] += len.min(hi) - lo;
                bin += 1;
            }
        }

        // the seeds are what mmseqs aligned and kept, read from the list
        // the run replayed
        let seeds = util::manifest::seeds_path(&root.join("results"), arm, &shard);
        let mut by_query: HashMap<String, Vec<String>> = HashMap::new();
        for line in BufReader::new(
            File::open(&seeds).with_context(|| format!("failed to open {}", seeds.display()))?,
        )
        .lines()
        {
            let line = line?;
            let mut f = line.split('\t');
            if let (Some(q), Some(t)) = (f.next(), f.next()) {
                by_query
                    .entry(q.to_string())
                    .or_default()
                    .push(t.to_string());
            }
        }

        pairs.sort_by(|a, b| a.query.cmp(&b.query));
        let mut queries: Vec<&str> = pairs.iter().map(|p| p.query.as_str()).collect();
        queries.extend(by_query.keys().map(String::as_str));
        queries.sort_unstable();
        queries.dedup();

        let mut at = 0;
        for query in queries {
            let ranks = match pf.query_key.get(query) {
                Some(&key) => pf.ranks(key)?,
                None => bail!("query {query:?} is not in {}", dir.display()),
            };
            let rank = |target: &str| -> Option<u64> {
                pf.target_key
                    .get(target)
                    .and_then(|k| ranks.get(k))
                    .map(|&r| r as u64)
            };

            for target in by_query.get(query).map(Vec::as_slice).unwrap_or_default() {
                match rank(target) {
                    Some(r) => Tally::bump(&mut tally.seeds, bin_of(r, width)),
                    None => bail!("seed {query} {target} is not in its prefilter list"),
                }
            }

            while at < pairs.len() && pairs[at].query == query {
                let p = &pairs[at];
                at += 1;
                match rank(&p.target) {
                    Some(r) => {
                        let bin = bin_of(r, width);
                        if p.kept {
                            Tally::bump(&mut tally.kept, bin);
                        }
                        if p.ceiling {
                            Tally::bump(&mut tally.ceiling, bin);
                        }
                        if p.kept && p.ceiling {
                            Tally::bump(&mut tally.hits, bin);
                            hits_table.row([
                                toil::Cell::from(shard.as_str()),
                                runs[run].name.as_str().into(),
                                p.query.as_str().into(),
                                p.target.as_str().into(),
                                r.into(),
                            ]);
                        }
                    }
                    None if p.kept => {
                        bail!(
                            "{} {} was kept without a prefilter entry",
                            p.query,
                            p.target
                        )
                    }
                    None => tally.ceiling_beyond += 1,
                }
            }
        }

        for bin in 0..tally.bins() {
            let (lo, hi) = range(bin, width);
            let get = |v: &Vec<u64>| v.get(bin).copied().unwrap_or(0);
            let (prefilter, seeds, kept, ceiling, hits) = (
                get(&tally.prefilter),
                get(&tally.seeds),
                get(&tally.kept),
                get(&tally.ceiling),
                get(&tally.hits),
            );
            table.row([
                toil::Cell::from(shard.as_str()),
                runs[run].name.as_str().into(),
                lo.into(),
                hi.into(),
                prefilter.into(),
                seeds.into(),
                kept.into(),
                ceiling.into(),
                hits.into(),
                frac(hits, prefilter).into(),
            ]);
        }
        table.row([
            toil::Cell::from(shard.as_str()),
            runs[run].name.as_str().into(),
            "-".into(),
            "-".into(),
            0u64.into(),
            0u64.into(),
            0u64.into(),
            tally.ceiling_beyond.into(),
            0u64.into(),
            "-".into(),
        ]);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))?;
    hits_table
        .write(hits_out)
        .with_context(|| format!("failed to write {}", hits_out.display()))?;
    lists_table
        .write(lists_out)
        .with_context(|| format!("failed to write {}", lists_out.display()))
}

fn frac(n: u64, of: u64) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_double_from_the_first_width() {
        assert_eq!(range(0, 200), (1, 200));
        assert_eq!(range(1, 200), (201, 400));
        assert_eq!(range(2, 200), (401, 800));
        for rank in [1, 200, 201, 400, 401, 800, 801] {
            let bin = bin_of(rank, 200);
            let (lo, hi) = range(bin, 200);
            assert!(lo <= rank && rank <= hi, "{rank} in bin {bin} {lo}..{hi}");
        }
    }
}
