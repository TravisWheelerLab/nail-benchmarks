//! How far down the prefilter list the hits keep coming.
//!
//! An arm seeded static and unbounded aligned everything its prefilter
//! returned, so every pair nail reported has a depth: its rank in its query's
//! prefilter list, which mmseqs writes in score order. Depth is binned in
//! ranges that double, the way prog's `n_take` does, and each bin counts the
//! prefilter pairs at that depth against the seeds, the pairs nail kept over
//! the cutoff, the pairs hmmer found, and the hits: kept and found both. The
//! hits over the prefilter pairs is what a stopping rule is betting on.
//!
//! The prefilter list is read out of the databases nail left under
//! `results/prefilter.<arm>.<shard>/`: mmseqs' index gives each query's
//! entry, split across one data file per thread whose offsets run on as if
//! the files were one, and nail's header databases give which key is which
//! name, in the order it wrote them.
//!
//! A family short of a cutoff for some tool is left out of every count and
//! every row, its prefilter list included. `#= left_out` on each table counts
//! the families and the pairs each run reported in them, and a run's last row
//! of `depth.tbl` carries its own count for the unit in `left_out`.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};

use crate::scores::runs::Reader as Runs;
use crate::scores::{LeftOut, Named};

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

    /// Every query in `keys` read once, on `threads`: how many lines its
    /// list holds, and the rank, from 1, of each target `wanted` asks about
    /// for it that the list holds.
    ///
    /// One pass rather than a map of the whole list: a list at -s 12 runs to
    /// millions of lines, and what is asked about it is a few thousand
    /// targets at most. A query `wanted` does not name is only counted.
    fn scan(
        &self,
        keys: &[u32],
        wanted: &HashMap<u32, HashSet<u32>>,
        threads: usize,
    ) -> anyhow::Result<Vec<(u32, Found)>> {
        use rayon::prelude::*;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .context("failed to build a thread pool")?;

        pool.install(|| {
            keys.par_iter()
                .map(|&key| {
                    let entry = self.entry(key)?;
                    let want = wanted.get(&key);
                    let mut found = Found {
                        len: memchr::memchr_iter(b'\n', &entry).count() as u64,
                        ranks: Vec::new(),
                    };

                    if let Some(want) = want {
                        // rank i + 1 for the i-th line, counting the ones
                        // after the data too, as the map of the whole list
                        // did: they come last, so no data line moves
                        for (i, line) in entry.split(|&b| b == b'\n').enumerate() {
                            if let Some(target) = leading_number(line)
                                && want.contains(&target)
                            {
                                found.ranks.push((target, i as u32 + 1));
                            }
                        }
                    }

                    Ok((key, found))
                })
                .collect()
        })
    }
}

/// What one scan of a query's list found.
struct Found {
    /// How many lines the list holds.
    len: u64,
    /// Each target asked about that the list holds, and its rank there.
    ranks: Vec<(u32, u32)>,
}

/// The number a prefilter line opens with, which is its target key, or
/// `None` for a line that opens with anything else: the `\0` that ends an
/// entry, or nothing.
fn leading_number(line: &[u8]) -> Option<u32> {
    let digits = line.iter().take_while(|b| b.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    let mut n: u32 = 0;
    for &b in &line[..digits] {
        n = n.checked_mul(10)?.checked_add((b - b'0') as u32)?;
    }
    Some(n)
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

/// How `pairs`, which one arm's seed list for one shard does not hold, were
/// kept out of it: the ones the prefilter never returned, and the ones it
/// returned past the rank where prog seeding stopped for their query, which
/// mmseqs was never asked to align. The rest it aligned and dropped.
///
/// prog seeding leaves each round's slice of the prefilter list under
/// `prog-seed/<round>/pdb`, the next `n_take` targets of every query still
/// active, so the pairs mmseqs aligned are the union of the slices. An arm
/// seeded static has no slices and aligned the whole list.
pub fn unseeded(
    root: &Path,
    arm: &str,
    shard: &str,
    pairs: &mut [(String, String)],
    threads: usize,
) -> anyhow::Result<(Vec<(String, String)>, Vec<(String, String)>)> {
    let dir = prefilter_dir(root, arm, shard);
    let pf = Prefilter::open(&dir)
        .with_context(|| format!("failed to open the prefilter database in {}", dir.display()))?;
    let rounds = slices(&dir.join("prog-seed"))?;

    // by query, so each query's entry is read once for all of its pairs
    pairs.sort_unstable();

    // which targets to look for in each query's list; a target the
    // database never keyed is missing before any list is read
    let mut wanted: HashMap<u32, HashSet<u32>> = HashMap::new();
    for (query, target) in pairs.iter() {
        let Some(&key) = pf.query_key.get(query) else {
            bail!("query {query:?} is not in {}", dir.display());
        };
        if let Some(&t) = pf.target_key.get(target) {
            wanted.entry(key).or_default().insert(t);
        }
    }
    let mut keys: Vec<u32> = wanted.keys().copied().collect();
    keys.sort_unstable();
    let listed: HashMap<u32, HashSet<u32>> = pf
        .scan(&keys, &wanted, threads)?
        .into_iter()
        .map(|(key, found)| (key, found.ranks.into_iter().map(|(t, _)| t).collect()))
        .collect();

    let (mut missing, mut stopped) = (Vec::new(), Vec::new());
    let mut at = 0;
    while at < pairs.len() {
        let query = pairs[at].0.as_str();
        let key = pf.query_key[query];
        let aligned = match rounds.is_empty() {
            true => None,
            false => Some(sliced(&rounds, key)?),
        };

        while at < pairs.len() && pairs[at].0 == query {
            let target = pf.target_key.get(&pairs[at].1);
            match target {
                Some(k) if listed.get(&key).is_some_and(|l| l.contains(k)) => {
                    if aligned.as_ref().is_some_and(|a| !a.contains(k)) {
                        stopped.push(pairs[at].clone());
                    }
                }
                _ => missing.push(pairs[at].clone()),
            }
            at += 1;
        }
    }

    Ok((missing, stopped))
}

/// One round's slice of the prefilter list: its index, and its one data file.
struct Slice {
    entries: HashMap<u32, (u64, u64)>,
    file: File,
}

/// Every round's slice under `prog-seed`, in round order; none where the
/// arm was not seeded prog.
fn slices(dir: &Path) -> anyhow::Result<Vec<Slice>> {
    if !dir.is_dir() {
        return Ok(Vec::new());
    }
    let mut rounds: Vec<(u32, PathBuf)> = std::fs::read_dir(dir)
        .with_context(|| format!("failed to read {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter_map(|p| Some((p.file_name()?.to_str()?.parse().ok()?, p)))
        .collect();
    rounds.sort_unstable();

    let mut out = Vec::with_capacity(rounds.len());
    for (_, round) in rounds {
        let pdb = round.join("pdb");
        out.push(Slice {
            entries: index(&round.join("pdb.index"))?
                .into_iter()
                .map(|(key, off, len)| (key, (off, len)))
                .collect(),
            file: File::open(&pdb).with_context(|| format!("failed to open {}", pdb.display()))?,
        });
    }
    Ok(out)
}

/// The target keys mmseqs was asked to align for one query, over every round.
fn sliced(rounds: &[Slice], key: u32) -> anyhow::Result<HashSet<u32>> {
    let mut keys = HashSet::new();
    for round in rounds {
        let Some(&(off, len)) = round.entries.get(&key) else {
            continue;
        };
        let mut entry = vec![0u8; len as usize];
        round.file.read_exact_at(&mut entry, off)?;
        for line in entry.split(|&b| b == b'\n') {
            let Some(field) = line.split(|&b| b == b'\t').next() else {
                continue;
            };
            if field.is_empty() || field == b"\0" {
                continue;
            }
            keys.insert(
                std::str::from_utf8(field)?.parse().with_context(|| {
                    format!("bad target key in a prog-seed slice for query {key}")
                })?,
            );
        }
    }
    Ok(keys)
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
    truth: Vec<u64>,
    /// Kept, and in hmmer's truth: the pairs a stopping rule is betting on.
    hits: Vec<u64>,
    /// hmmer's pairs the prefilter never returned, so they have no depth.
    truth_beyond: u64,
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
            &self.truth,
            &self.hits,
        ]
        .iter()
        .map(|v| v.len())
        .max()
        .unwrap_or(0)
    }
}

/// A pair one run kept, or hmmer found, waiting on its rank.
struct Pair {
    query: String,
    target: String,
    kept: bool,
    truth: bool,
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
    threads: usize,
) -> anyhow::Result<()> {
    ensure!(width > 0, "--bin must be at least 1");

    let mut scores = Runs::open(table)?;
    let hmmer = scores.meta().hmmer()?;
    let runs = scores.meta().runs.clone();
    let named = Named::read(&scores.meta().cutoffs, scores.meta().c)?;

    let mut left = LeftOut::new(runs.len());
    let mut left_of: HashMap<(usize, String), u64> = HashMap::new();

    // every pair that hmmer found or a run kept, grouped by the arm and the
    // shard whose prefilter list says where it sat
    let mut pairs: HashMap<(usize, String), Vec<Pair>> = HashMap::new();
    // what hmmer found per unit, the denominator a sensitivity is read over
    let mut truth_of: indexmap::IndexMap<String, u64> = indexmap::IndexMap::new();

    scores.each(|row| {
        if !row.row().complete() {
            for run in 0..runs.len() {
                *left_of
                    .entry((run, row.row().shard().to_string()))
                    .or_default() += u64::from(row.present(run));
            }
            let family = std::str::from_utf8(row.row().field(0))
                .context("a query name that is not utf-8")?;
            left.pair(family, |run| row.present(run));
            return Ok(());
        }

        let truth = row.hmmer_found(hmmer);
        if truth {
            *truth_of.entry(row.row().shard().to_string()).or_default() += 1;
        }
        let query = String::from_utf8_lossy(row.row().field(0)).into_owned();
        let target = String::from_utf8_lossy(row.row().field(1)).into_owned();

        for (run, column) in runs.iter().enumerate() {
            if run == hmmer || column.seeds.is_none() {
                continue;
            }
            let kept = row.present(run) && row.row().passed(run);
            if !kept && !truth {
                continue;
            }
            pairs
                .entry((run, row.row().shard().to_string()))
                .or_default()
                .push(Pair {
                    query: query.clone(),
                    target: target.clone(),
                    kept,
                    truth,
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
        toil::Column::new("truth"),
        toil::Column::new("hits"),
        toil::Column::new("frac").fixed(4),
        toil::Column::new("left_out"),
    ]));
    table.meta("bin", [width.to_string()]);
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
    for (unit, n) in &truth_of {
        table.meta("truth", [unit.clone(), n.to_string()]);
        hits_table.meta("truth", [unit.clone(), n.to_string()]);
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
        queries.extend(
            by_query
                .keys()
                .map(String::as_str)
                .filter(|query| named.complete(query)),
        );
        queries.sort_unstable();
        queries.dedup();

        // every list is read once, for its length, and the targets the
        // seeds and the pairs name are looked up on the way through
        let mut wanted: HashMap<u32, HashSet<u32>> = HashMap::new();
        for query in &queries {
            let Some(&key) = pf.query_key.get(*query) else {
                bail!("query {query:?} is not in {}", dir.display());
            };
            let want = wanted.entry(key).or_default();
            for target in by_query.get(*query).map(Vec::as_slice).unwrap_or_default() {
                want.extend(pf.target_key.get(target));
            }
        }
        for p in &pairs {
            if let (Some(&key), Some(&t)) =
                (pf.query_key.get(&p.query), pf.target_key.get(&p.target))
            {
                wanted.entry(key).or_default().insert(t);
            }
        }

        // in the order the data holds the lists, which is how they are
        // read fastest and how the per-query table has always run
        let mut keys: Vec<(u64, u64, u32)> = pf
            .entries
            .iter()
            .map(|(&key, &(off, len))| (off, len, key))
            .collect();
        keys.sort_unstable();

        let name_of: HashMap<u32, &str> = pf
            .query_key
            .iter()
            .map(|(name, &key)| (key, name.as_str()))
            .collect();

        // a key with no name cannot be judged, and is kept as it always was
        let keys: Vec<u32> = keys
            .into_iter()
            .map(|(_, _, key)| key)
            .filter(|key| match name_of.get(key) {
                Some(name) if !named.complete(name) => {
                    left.family(name);
                    false
                }
                _ => true,
            })
            .collect();

        let mut found: HashMap<u32, HashMap<u32, u32>> = HashMap::new();
        for (key, list) in pf.scan(&keys, &wanted, threads)? {
            // the prefilter pairs per bin come from the length of every
            // query's list: a list of L holds min(L, hi) - lo + 1 of each bin
            lists_table.row([
                toil::Cell::from(shard.as_str()),
                runs[run].name.as_str().into(),
                (*name_of.get(&key).unwrap_or(&"?")).into(),
                list.len.into(),
            ]);
            let mut bin = 0;
            loop {
                let (lo, hi) = range(bin, width);
                if list.len < lo {
                    break;
                }
                Tally::bump(&mut tally.prefilter, bin);
                tally.prefilter[bin] += list.len.min(hi) - lo;
                bin += 1;
            }
            found.insert(key, list.ranks.into_iter().collect());
        }

        let mut at = 0;
        for query in queries {
            let ranks = &found[&pf.query_key[query]];
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
                        if p.truth {
                            Tally::bump(&mut tally.truth, bin);
                        }
                        if p.kept && p.truth {
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
                    None => tally.truth_beyond += 1,
                }
            }
        }

        for bin in 0..tally.bins() {
            let (lo, hi) = range(bin, width);
            let get = |v: &Vec<u64>| v.get(bin).copied().unwrap_or(0);
            let (prefilter, seeds, kept, truth, hits) = (
                get(&tally.prefilter),
                get(&tally.seeds),
                get(&tally.kept),
                get(&tally.truth),
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
                truth.into(),
                hits.into(),
                frac(hits, prefilter).into(),
                "-".into(),
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
            tally.truth_beyond.into(),
            0u64.into(),
            "-".into(),
            left_of
                .get(&(run, shard.clone()))
                .copied()
                .unwrap_or(0)
                .into(),
        ]);
    }

    // after the groups, since the prefilter lists add families the rows
    // never named
    for words in left.lines(&runs) {
        table.meta("left_out", words.clone());
        hits_table.meta("left_out", words.clone());
        lists_table.meta("left_out", words);
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
