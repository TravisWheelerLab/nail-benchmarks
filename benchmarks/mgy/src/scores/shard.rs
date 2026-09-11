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
use libsail::tbl::hmmer::{HmmerDomTable, HmmerTable};
use libsail::tbl::nail::NailTable;

use util::manifest;

use super::scan::{self, Lines};
use super::{Column, Cutoffs, Queries, Tool};

/// How many bits of a key the target holds.
const TID: u32 = 41;

/// One pair as one run reported it.
struct Hit {
    key: u64,
    score: f32,
    run: u16,
    /// hmmer's inclusion count. Zero for every other tool, which do not break
    /// a hit into domains and so have nothing to include.
    inc: u16,
}

/// One domain of one hmmer hit, in the order the domtbl listed it.
struct Dom {
    key: u64,
    score: f32,
    ord: u32,
}

/// One pair, and everything a row could say about it.
pub struct Pair<'a> {
    pub query: u32,
    pub target: &'a [u8],
    /// One character per run, in ledger order: the tool's letter, uppercase
    /// where that run reported the pair at or above its family's cutoff.
    pub pass: &'a [u8],
    /// One per tool, at [`Tool::at`], `None` where no run of that tool
    /// reported the pair.
    pub scores: &'a [Option<f32>; 3],
    /// hmmer's inclusion count, `None` where hmmer did not report the pair.
    pub inc: Option<u16>,
    /// hmmer's domain scores, in domtbl order.
    pub doms: &'a [f32],
}

/// What one shard came to.
#[derive(Clone, Copy, Debug, Default)]
pub struct Count {
    pub rows: u64,
    pub hits: u64,
    /// Pairs two runs of one tool gave different scores, which the table has
    /// one column for and so cannot show both of.
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
    doms: Vec<Dom>,
    pass: Vec<u8>,
    dom_scores: Vec<f32>,
    name: Vec<u8>,
}

/// Where a shard's results are, and what to make of them.
pub struct Shard<'a> {
    pub results: &'a Path,
    pub runs: &'a [Column],
    pub queries: &'a Queries,
    pub cutoffs: &'a Cutoffs,
    /// Which run's `.domtbl` carries the domain breakdown.
    pub hmmer: Option<usize>,
}

impl Shard<'_> {
    /// Read every run's table for one shard and hand each pair to `render`.
    ///
    /// A pair earns a row by clearing some run's cutoff or by hmmer having
    /// reported it; the rest are read, folded and dropped.
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
        scratch.doms.clear();

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
                Tool::Hmmer => rows.hmmer(&path)?,
            };

            if !read {
                return Ok(None);
            }

            if self.hmmer == Some(at) {
                let path = manifest::dom_path(self.results, &column.run.name, shard);
                if !doms(&path, &mut keys, self.queries, &mut scratch.doms)? {
                    return Ok(None);
                }
            }
        }

        if let Keys::Named(names) = &mut keys {
            // interned in the order met, so the keys have to be put back into
            // the order the names sort before anything is sorted by them
            let rank = names.order();
            for hit in &mut scratch.hits {
                hit.key = rerank(hit.key, &rank);
            }
            for dom in &mut scratch.doms {
                dom.key = rerank(dom.key, &rank);
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
            doms,
            pass,
            dom_scores,
            name,
        } = scratch;

        hits.sort_unstable_by_key(|hit| (hit.key, hit.run));
        doms.sort_unstable_by_key(|dom| (dom.key, dom.ord));

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
        let mut dom_at = 0usize;

        while at < hits.len() {
            let key = hits[at].key;
            let query = (key >> TID) as u32;

            pass.clear();
            pass.extend_from_slice(&letters);

            let mut scores: [Option<f32>; 3] = [None; 3];
            let mut inc: Option<u16> = None;

            while at < hits.len() && hits[at].key == key {
                let run = hits[at].run;

                // a run can report a pair more than once; the best of them is
                // the one a threshold would see
                let mut best = f32::NEG_INFINITY;
                let mut included = 0u16;

                while at < hits.len() && hits[at].key == key && hits[at].run == run {
                    best = best.max(hits[at].score);
                    included = included.max(hits[at].inc);
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

                match scores[tool.at()] {
                    None => scores[tool.at()] = Some(best),
                    // runs of one tool agree on a pair's score by
                    // construction, so a disagreement is worth counting rather
                    // than picking between
                    Some(first) if first != best => count.disagreements += 1,
                    Some(_) => {}
                }

                if tool == Tool::Hmmer {
                    inc = Some(included);
                }
            }

            // every domain of this pair, in the order the domtbl listed them
            dom_scores.clear();
            while dom_at < doms.len() && doms[dom_at].key < key {
                dom_at += 1;
            }
            while dom_at < doms.len() && doms[dom_at].key == key {
                dom_scores.push(doms[dom_at].score);
                dom_at += 1;
            }

            let passed = pass.iter().any(u8::is_ascii_uppercase);
            if !passed && scores[Tool::Hmmer.at()].is_none() {
                continue;
            }

            keys.name(key, name);
            render(&Pair {
                query,
                target: &name[..],
                pass: &pass[..],
                scores: &scores,
                inc,
                doms: &dom_scores[..],
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

        while let Some((at, line)) = lines.next()? {
            if !checked {
                if !scan::fits::<C>(line) {
                    bail!(
                        "{}:{at} has {} fields, this layout writes {}",
                        path.display(),
                        scan::count(line),
                        C::N_FIELDS
                    );
                }
                checked = true;
            }

            let Some([query, target, score]) = scan::hit::<C>(line) else {
                bail!("{}:{at} is short of fields", path.display());
            };

            if !self.push(query, target, score, 0, path, at)? {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// hmmer's `--tblout`, which carries one column the others do not.
    fn hmmer(&mut self, path: &Path) -> anyhow::Result<bool> {
        let mut lines = open(path)?;
        let mut checked = false;

        while let Some((at, line)) = lines.next()? {
            if !checked {
                if !scan::fits::<HmmerTable>(line) {
                    bail!(
                        "{}:{at} has {} fields, hmmer's --tblout writes {}",
                        path.display(),
                        scan::count(line),
                        HmmerTable::N_FIELDS
                    );
                }
                checked = true;
            }

            let Some([query, target, score, inc]) = scan::hmmer_hit(line) else {
                bail!("{}:{at} is short of fields", path.display());
            };

            let count: u16 = std::str::from_utf8(inc)
                .ok()
                .and_then(|text| text.parse().ok())
                .with_context(|| {
                    format!(
                        "{}:{at} has an inc column of {:?}",
                        path.display(),
                        String::from_utf8_lossy(inc)
                    )
                })?;

            if !self.push(query, target, score, count, path, at)? {
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
        inc: u16,
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
                    format!("{}:{} has a query name that is not text", path.display(), line)
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
            inc,
        });

        Ok(true)
    }
}

/// Every domain of every hit in one `--domtblout`, in the order it listed them.
fn doms(
    path: &Path,
    keys: &mut Keys,
    queries: &Queries,
    out: &mut Vec<Dom>,
) -> anyhow::Result<bool> {
    let mut lines = open(path)?;
    let mut checked = false;
    let mut last: Option<(Vec<u8>, u32)> = None;
    let mut ord = 0u32;

    while let Some((at, line)) = lines.next()? {
        if !checked {
            if !scan::fits::<HmmerDomTable>(line) {
                bail!(
                    "{}:{at} has {} fields, hmmer's --domtblout writes {}",
                    path.display(),
                    scan::count(line),
                    HmmerDomTable::N_FIELDS
                );
            }
            checked = true;
        }

        let Some([query, target, score]) = scan::hit::<HmmerDomTable>(line) else {
            bail!("{}:{at} is short of fields", path.display());
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
                        "{}:{at} reports family {name:?}, which is not in the query set",
                        path.display()
                    )
                })?;

                last = Some((query.to_vec(), id));
                id
            }
        };

        let score = scan::score(score).with_context(|| {
            format!(
                "{}:{at} has a domain score of {:?}",
                path.display(),
                String::from_utf8_lossy(score)
            )
        })?;

        out.push(Dom {
            key: (qid as u64) << TID | tid,
            score,
            ord,
        });

        ord += 1;
    }

    Ok(true)
}

fn open(path: &Path) -> anyhow::Result<Lines<File>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    Ok(Lines::new(file))
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
