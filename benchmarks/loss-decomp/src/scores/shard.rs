//! Every run's hits for one shard, folded into one row per pair.
//!
//! This is the part no analysis owns. A shard's tables are read into one flat
//! vector of hits, sorted once, and walked in a single pass that hands each
//! pair to whatever is writing the table. What the row looks like belongs to
//! the caller; what a pair is does not.
//!
//! A flat vector rather than a map: the work is bandwidth, one structure to
//! fill and sort beats millions of allocations, and a target name that is an
//! `MGYP` number is a key already -- so nothing here hashes a string per row.
//!
//! ```text
//! key = qid << 41 | tid
//! ```
//!
//! `tid` is the number in the target's name, under 2^40, and `qid` the
//! family's rank among the query names. Both order as their names do, so
//! sorting the key sorts the pairs.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use anyhow::{Context, bail};

use libsail::tbl::HitColumns;
use libsail::tbl::blast::BlastTable;
use libsail::tbl::nail::NailTable;

use util::manifest;

use crate::scores::scan;
use crate::scores::{Column, Cutoffs, Queries, Tool};
use libsail::lines::Rows as TableRows;

/// How many bits of a key the target holds.
const TID: u32 = 41;

/// One pair as one run reported it.
struct Hit {
    key: u64,
    score: f32,
    run: u16,
}

/// One pair, and everything a row could say about it.
pub struct Pair<'a> {
    pub query: u32,
    pub target: &'a [u8],
    /// One character per run, in ledger order: the tool's letter, uppercase
    /// where that run reported the pair at or above its family's cutoff, and
    /// `-` where the pair was not in that run's seed list at all. Not seeded
    /// and seeded but unreported are different answers.
    pub pass: &'a [u8],
    /// One per run, in ledger order, `None` where that run did not report
    /// the pair.
    pub scores: &'a [Option<f32>],
}

/// What one shard came to.
#[derive(Clone, Copy, Debug, Default)]
pub struct Count {
    pub rows: u64,
    pub hits: u64,
    /// Pairs two runs of one tool gave different scores, for a table that
    /// folds them into one column.
    pub disagreements: u64,
}

impl Count {
    pub fn add(&mut self, other: Count) {
        self.rows += other.rows;
        self.hits += other.hits;
        self.disagreements += other.disagreements;
    }
}

/// The buffers a shard is read into, kept across shards so a worker allocates
/// once rather than per shard.
#[derive(Default)]
pub struct Scratch {
    hits: Vec<Hit>,
    pass: Vec<u8>,
    scores: Vec<Option<f32>>,
    seeds: Vec<Vec<u64>>,
    name: Vec<u8>,
}

/// Where a shard's results are, and what to make of them.
pub struct Shard<'a> {
    pub results: &'a Path,
    pub runs: &'a [Column],
    pub queries: &'a Queries,
    pub cutoffs: &'a Cutoffs,
    /// Which seed list each run replayed, as an index into [`Self::lists`],
    /// and `None` for a run that replayed none.
    ///
    /// Many runs to one list: a pruning sweep replays a single seeding into
    /// every cell, a seeding sweep gives every arm its own. The seed list is
    /// read for a character on rows, not for rows of its own -- a pair nothing
    /// reported is not among the hits and never reaches a renderer.
    pub seeds: &'a [Option<usize>],
    /// The distinct seed lists, by the name a run's `seeds` setting gives.
    pub lists: &'a [String],
}

impl Shard<'_> {
    /// Read every run's table for one shard and hand each pair to `render`.
    ///
    /// A pair earns a row by clearing some run's cutoff; the rest are read,
    /// folded and dropped.
    pub fn collect(
        &self,
        shard: &str,
        scratch: &mut Scratch,
        render: &mut impl FnMut(&Pair<'_>) -> anyhow::Result<()>,
    ) -> anyhow::Result<Count> {
        // a shard whose targets are not all MGnify names is read again with
        // its names interned, since by the time one is met the ones before it
        // are numbers and their names are gone
        match self.pass(shard, Keys::Mgyp, scratch, render)? {
            Some(count) => Ok(count),
            None => self
                .pass(shard, Keys::Named(Names::default()), scratch, render)?
                .context("a shard of interned names still would not key"),
        }
    }

    /// One attempt at a shard. `None` where the keys could not hold a name.
    fn pass(
        &self,
        shard: &str,
        mut keys: Keys,
        scratch: &mut Scratch,
        render: &mut impl FnMut(&Pair<'_>) -> anyhow::Result<()>,
    ) -> anyhow::Result<Option<Count>> {
        scratch.hits.clear();
        scratch.seeds.clear();

        // before the hit tables rather than after: the seed list holds every
        // pair any nail run could report, so a shard whose names will not key
        // as MGYP is detected here, and the retry re-reads one small file
        // than every table in the shard
        scratch.seeds.resize(self.lists.len(), Vec::new());
        for (list, into) in self.lists.iter().zip(scratch.seeds.iter_mut()) {
            into.clear();
            let path = manifest::seeds_path(self.results, list, shard);
            if !seeds(&path, &mut keys, self.queries, into)? {
                return Ok(None);
            }
        }

        for (at, column) in self.runs.iter().enumerate() {
            if !column.shards.iter().any(|covered| covered == shard) {
                continue;
            }

            let path = manifest::table_path(self.results, &column.run.name, shard);
            let mut rows = Rows {
                keys: &mut keys,
                queries: self.queries,
                hits: &mut scratch.hits,
                last: None,
                at: at as u16,
            };

            let read = match column.run.tool {
                Tool::Nail => rows.table::<NailTable>(&path)?,
                Tool::Mmseqs => rows.table::<BlastTable>(&path)?,
            };

            if !read {
                return Ok(None);
            }
        }

        if let Keys::Named(names) = &mut keys {
            // interned in the order met, so the keys have to be put back into
            // the order the names sort before anything is sorted by them
            let rank = names.order();
            for hit in &mut scratch.hits {
                hit.key = rerank(hit.key, &rank);
            }
            for list in &mut scratch.seeds {
                for key in list.iter_mut() {
                    *key = rerank(*key, &rank);
                }
            }
        }

        Ok(Some(self.fold(&keys, scratch, render)?))
    }

    /// Walk the sorted hits once, a pair at a time.
    fn fold(
        &self,
        keys: &Keys,
        scratch: &mut Scratch,
        render: &mut impl FnMut(&Pair<'_>) -> anyhow::Result<()>,
    ) -> anyhow::Result<Count> {
        let Scratch {
            hits,
            pass,
            scores,
            seeds,
            name,
        } = scratch;

        hits.sort_unstable_by_key(|hit| (hit.key, hit.run));
        for list in seeds.iter_mut() {
            list.sort_unstable();
        }

        let mut count = Count {
            hits: hits.len() as u64,
            ..Count::default()
        };

        let letters: Vec<u8> = self
            .runs
            .iter()
            .map(|column| column.run.tool.letter())
            .collect();

        let mut at = 0usize;

        // one cursor per list rather than per run: several runs replay one
        // seeding, and advancing per run would walk a shared list twice
        let mut seed_at = vec![0usize; seeds.len()];
        let mut holds = vec![false; seeds.len()];

        while at < hits.len() {
            let key = hits[at].key;
            let query = (key >> TID) as u32;

            pass.clear();
            pass.extend_from_slice(&letters);

            // cleared rather than resized: `resize` on a vector that already
            // has the length leaves the previous pair's scores in it, which
            // would invent one for a run that reported nothing
            scores.clear();
            scores.resize(self.runs.len(), None);

            while at < hits.len() && hits[at].key == key {
                let run = hits[at].run;

                // a run can report a pair more than once; the best of them is
                // the one a threshold would see
                let mut best = f32::NEG_INFINITY;

                while at < hits.len() && hits[at].key == key && hits[at].run == run {
                    best = best.max(hits[at].score);
                    at += 1;
                }

                let tool = self.runs[run as usize].run.tool;

                if self
                    .cutoffs
                    .get(tool, query)
                    .is_some_and(|cutoff| best >= cutoff)
                {
                    pass[run as usize] = pass[run as usize].to_ascii_uppercase();
                }

                scores[run as usize] = Some(best);
            }

            // the lists are sorted by the same key the hits are, so each is
            // walked forward once across the whole shard rather than searched
            // per pair. a seed list with the pair twice is harmless
            for (li, list) in seeds.iter().enumerate() {
                while seed_at[li] < list.len() && list[seed_at[li]] < key {
                    seed_at[li] += 1;
                }
                holds[li] = list.get(seed_at[li]) == Some(&key);
            }

            // a run whose seeding never offered the pair did not miss it, and
            // the character says so rather than a column shared by every run
            for (run, list) in self.seeds.iter().enumerate() {
                if list.is_some_and(|li| !holds[li]) {
                    pass[run] = b'-';
                }
            }

            if !pass.iter().any(u8::is_ascii_uppercase) {
                continue;
            }

            keys.name(key, name);
            render(&Pair {
                query,
                target: &name[..],
                pass: &pass[..],
                scores: &scores[..],
            })?;

            count.rows += 1;
        }

        Ok(count)
    }
}

/// One results table's rows, as they go into the shard's hits.
struct Rows<'a> {
    keys: &'a mut Keys,
    queries: &'a Queries,
    hits: &'a mut Vec<Hit>,
    /// The last query name and its id. A table's rows come in runs of one
    /// query, so this answers most of the lookups without a hash.
    last: Option<(Vec<u8>, u32)>,
    at: u16,
}

impl Rows<'_> {
    /// A table in layout `C`. `false` where a target name did not key.
    fn table<C: HitColumns>(&mut self, path: &Path) -> anyhow::Result<bool> {
        let mut lines = open(path)?;
        let mut checked = false;

        while lines.advance()? {
            let (at, line) = (lines.line(), lines.row());
            if !checked {
                if !scan::fits::<C>(line) {
                    bail!(
                        "{}:{at} has {} fields, this layout writes {}",
                        path.display(),
                        scan::count(line),
                        C::N_COLUMNS
                    );
                }
                checked = true;
            }

            let Some([query, target, score]) = scan::hit::<C>(line) else {
                bail!("{}:{at} is short of fields", path.display());
            };

            if !self.push(query, target, score, path, at)? {
                return Ok(false);
            }
        }

        Ok(true)
    }

    fn push(
        &mut self,
        query: &[u8],
        target: &[u8],
        score: &[u8],
        path: &Path,
        line: u64,
    ) -> anyhow::Result<bool> {
        let Some(tid) = self.keys.tid(target) else {
            return Ok(false);
        };

        let qid = match &self.last {
            Some((name, id)) if name == query => *id,
            _ => {
                let name = std::str::from_utf8(query).with_context(|| {
                    format!(
                        "{}:{} has a query name that is not text",
                        path.display(),
                        line
                    )
                })?;

                let id = self.queries.id(name).with_context(|| {
                    format!(
                        "{}:{} reports family {name:?}, which is not in the query set",
                        path.display(),
                        line
                    )
                })?;

                self.last = Some((query.to_vec(), id));
                id
            }
        };

        let score = scan::score(score).with_context(|| {
            format!(
                "{}:{} has a score of {:?}",
                path.display(),
                line,
                String::from_utf8_lossy(score)
            )
        })?;

        self.hits.push(Hit {
            key: (qid as u64) << TID | tid,
            score,
            run: self.at,
        });

        Ok(true)
    }
}

/// The (query, target) pairs one shard's seeding found.
///
/// nail's `--seeds-out` is its own two-column list -- profile then sequence,
/// whitespace-separated, no header -- rather than a hit table, so it is read
/// here rather than through a `libsail` layout. Keyed the same way the hits
/// are, since it is walked beside them.
fn seeds(
    path: &Path,
    keys: &mut Keys,
    queries: &Queries,
    out: &mut Vec<u64>,
) -> anyhow::Result<bool> {
    let mut lines = open(path)?;
    let mut last: Option<(Vec<u8>, u32)> = None;

    while lines.advance()? {
        let (at, line) = (lines.line(), lines.row());
        let Some([query, target]) = scan::fields(line, [0, 1]) else {
            bail!(
                "{}:{at} has {} fields, a seed list is a query and a target",
                path.display(),
                scan::count(line)
            );
        };

        let Some(tid) = keys.tid(target) else {
            return Ok(false);
        };

        let qid = match &last {
            Some((name, id)) if name == query => *id,
            _ => {
                let name = std::str::from_utf8(query).with_context(|| {
                    format!("{}:{at} has a query name that is not text", path.display())
                })?;

                let id = queries.id(name).with_context(|| {
                    format!(
                        "{}:{at} seeds family {name:?}, which is not in the query set",
                        path.display()
                    )
                })?;

                last = Some((query.to_vec(), id));
                id
            }
        };

        out.push((qid as u64) << TID | tid);
    }

    Ok(true)
}

fn open(path: &Path) -> anyhow::Result<TableRows<File>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    // a comment, and a line with nothing but spaces on it. Rows drops an empty
    // line on its own, and judges the rest on the first byte: no hit table
    // here opens a row with whitespace, so rejecting one is the blank-line
    // skip this reader has always had
    Ok(TableRows::new(file, |b| {
        b == b'#' || b.is_ascii_whitespace()
    }))
}

// --------------------------------------------------------------------- keys

/// How a shard's target names become the numbers a key is built from.
enum Keys {
    /// `MGYP` and twelve digits, which is a number already.
    Mgyp,
    /// Anything else, interned.
    Named(Names),
}

#[derive(Default)]
struct Names {
    at: HashMap<Vec<u8>, u32>,
    names: Vec<Vec<u8>>,
}

impl Names {
    /// Put the names in order, and say where each one moved to.
    ///
    /// Reading is over by the time this is called, so the map from name to id
    /// is left behind: what is wanted from here on is the name at a rank.
    fn order(&mut self) -> Vec<u64> {
        let mut order: Vec<u32> = (0..self.names.len() as u32).collect();
        order.sort_unstable_by(|&a, &b| self.names[a as usize].cmp(&self.names[b as usize]));

        let mut rank = vec![0u64; order.len()];
        let mut sorted = Vec::with_capacity(order.len());

        for (to, &from) in order.iter().enumerate() {
            rank[from as usize] = to as u64;
            sorted.push(std::mem::take(&mut self.names[from as usize]));
        }

        self.names = sorted;
        self.at = HashMap::new();

        rank
    }
}

impl Keys {
    /// The target half of a key, or `None` where this shard's names cannot go
    /// in the keys as they stand.
    fn tid(&mut self, target: &[u8]) -> Option<u64> {
        match self {
            Keys::Mgyp => scan::mgyp(target),
            Keys::Named(names) => Some(match names.at.get(target) {
                Some(&id) => id as u64,
                None => {
                    let id = names.names.len() as u64;
                    names.at.insert(target.to_vec(), id as u32);
                    names.names.push(target.to_vec());
                    id
                }
            }),
        }
    }

    /// The name a key's target half stands for, written into `out`.
    fn name(&self, key: u64, out: &mut Vec<u8>) {
        let tid = key & ((1 << TID) - 1);

        out.clear();
        match self {
            Keys::Mgyp => {
                out.extend_from_slice(b"MGYP");
                let mut digits = [b'0'; 12];
                let mut left = tid;
                for digit in digits.iter_mut().rev() {
                    *digit = b'0' + (left % 10) as u8;
                    left /= 10;
                }
                out.extend_from_slice(&digits);
            }
            // the keys were put into name order, and so were the names
            Keys::Named(names) => out.extend_from_slice(&names.names[tid as usize]),
        }
    }
}

/// One key's target half, moved from the order it was met to the order it
/// sorts.
fn rerank(key: u64, rank: &[u64]) -> u64 {
    let mask = (1u64 << TID) - 1;
    (key & !mask) | rank[(key & mask) as usize]
}
