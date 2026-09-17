//! `runs.tbl`: one row per pair, one score column per run.
//!
//! This is the table for the pipelines that move one tool's parameters and
//! ask what that costs. recall's `scores.tbl` keeps one column per tool on
//! the premise that a tool scores a pair the same wherever it reports it;
//! `-A` and `-B` are exactly what break that premise, since pruning
//! constrains the dynamic programming and a pair can come back with a lower
//! score. So the score is kept per run here, and the price -- a column per
//! cell rather than per tool -- is one this pipeline can pay: it searches one
//! shard, where recall searches a thousand.
//!
//! ```text
//! #= format runs 1
//! #= query <count> <residues> <bytes>
//! #= target <shard> <count> <residues> <bytes>
//! #= seed <seeding> <wall_s>
//! #= cutoffs <path> c=<n>
//! #= run <name> <tool> <wall_s> [k=v ...]
//! #= pass <run name> ...
//! # query target           pass seeded A2.0-B4.0 ... full   hmmer  inc dom
//! # ----- ---------------- ---- ------ --------- --- ------ ------ --- ---
//! #= shard 1
//! 2-Hacid_dh_C MGYP000522683479 NNH  y      98.8      ... 98.8   98.7   1   98.1
//! #= end <rows>
//! ```
//!
//! `seeded` is `y` or `n`, and the column is there exactly when the pipeline
//! kept a seed list beside its results. Seeding nothing and never having
//! seeded are different answers, so a pipeline that never seeded has no column
//! rather than a column of `n`.
//!
//! A `-` in a run's column means that run did not report the pair. Absent is
//! not a score: a zero or a NaN would compare against a threshold and look
//! like one, which for this table is the whole question.

use std::fmt::Write as _;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, ensure};

use tabl::{Column, Schema, Stream, Widths};

use util::ledger::{self, Ledger};
use util::set::Set;

use crate::collect::{self, Job};
use crate::frame::{Frame, Layout};
use crate::shard::{Count, Pair, Scratch, Shard};
use crate::{Cutoffs, Meta, Queries, label, runs};

/// The `#= format` line a runs table opens with.
pub const FORMAT: &str = "#= format runs 1";

/// How wide a target name is written. MGnify's are sixteen characters, and a
/// wider one overruns rather than widening the column, since the header has
/// gone out before the first row is read.
const TARGET: usize = 16;

/// How wide a score is written: four digits, a point and a place.
const SCORE: usize = 6;

/// Where the score columns start: after query, target and pass.
const SCORES_AT: usize = 3;

pub struct Args<'a> {
    /// The pipeline directory: `ledger.tbl` and `results/`.
    pub dir: &'a Path,
    pub query_hmm: &'a Path,
    /// The set that was searched, for what each unit came to.
    pub set: &'a Set,
    pub cutoffs: &'a Path,
    pub c: usize,
    pub out: &'a Path,
    pub threads: usize,
    /// How many bytes the collectors may hold between them.
    pub mem: u64,
}

/// Read a finished pipeline directory into `runs.tbl`.
pub fn collect(args: Args<'_>) -> anyhow::Result<Count> {
    let ran = Ledger::load(args.dir)?;
    ledger::warn(ran.failed(), "command(s)");

    let columns = runs(&ran)?;
    ensure!(
        !columns.is_empty(),
        "no finished runs in {}",
        args.dir.display()
    );

    let shards = ran.shards();
    ensure!(!shards.is_empty(), "no shards in {}", args.dir.display());

    // what the pipeline itself recorded about seeding: a stage row per seed
    // list. Keyed by the list rather than by the shard, since a seeding sweep
    // writes one per arm and the shard cannot tell them apart
    // summed over the shards a seeding covered, the way a run's wall clock is:
    // one line per seeding rather than one per seeding per shard
    let seeds: Vec<(String, f64)> = {
        let mut totals: Vec<(String, f64)> = Vec::new();
        for row in ran.stage("seed") {
            let name = row
                .params
                .get(util::manifest::SEEDS)
                .cloned()
                .unwrap_or_else(|| row.shard.clone());

            let wall = row.wall_s.unwrap_or_default();
            match totals.iter_mut().find(|(at, _)| *at == name) {
                Some((_, total)) => *total += wall,
                None => totals.push((name, wall)),
            }
        }
        totals
    };

    let queries = Queries::from_hmm(args.query_hmm)?;
    let cutoffs = Cutoffs::read(args.cutoffs, args.c, &queries)?;

    let meta = Meta {
        query: queries.size,
        targets: crate::target_sizes(args.set, &shards)?,
        seeds,
        cutoffs: collect::absolute(args.cutoffs),
        c: args.c,
        runs: columns.iter().map(|column| column.run.clone()).collect(),
        tools: ran.tools().to_vec(),
    };

    // every run here is measured against what hmmer found, and the domain
    // breakdown is hmmer's alone, so which run is hmmer's has to be a single
    // answer
    let hmmer = meta.hmmer()?;

    // one entry per run, naming the seeding it replayed, and the distinct
    // seedings behind them
    let lists: Vec<String> = {
        let mut seen: Vec<String> = Vec::new();
        for column in &columns {
            if let Some(name) = &column.run.seeds
                && !seen.contains(name)
            {
                seen.push(name.clone());
            }
        }
        seen
    };

    let seed_of: Vec<Option<usize>> = columns
        .iter()
        .map(|column| {
            column
                .run
                .seeds
                .as_ref()
                .and_then(|name| lists.iter().position(|list| list == name))
        })
        .collect();

    let schema = schema(&meta);
    let widths = schema.widths();

    let file = std::fs::File::create(args.out)
        .with_context(|| format!("failed to create {}", args.out.display()))?;
    let mut out = BufWriter::with_capacity(1 << 20, file);

    meta.write(FORMAT, &mut out)?;
    Stream::new(schema.clone(), widths.clone(), &mut out).header()?;

    let results = args.dir.join("results");
    let work = Shard {
        results: &results,
        runs: &columns,
        queries: &queries,
        cutoffs: &cutoffs,
        hmmer: Some(hmmer),
        seeds: &seed_of,
        lists: &lists,
    };

    let job = Job {
        work: &work,
        shards: &shards,
        threads: args.threads,
        mem: args.mem,
    };

    let count = collect::blocks(
        job,
        |work, shard, scratch| block(work, shard, scratch, &schema, &widths),
        &mut out,
    )?;

    writeln!(out, "#= end {}", count.rows)?;
    out.flush()?;

    Ok(count)
}

/// One shard's rows, rendered into memory.
fn block(
    work: &Shard<'_>,
    shard: &str,
    scratch: &mut Scratch,
    schema: &Schema,
    widths: &Widths,
) -> anyhow::Result<(Vec<u8>, Count)> {
    let mut out = Stream::continued(schema.clone(), widths.clone(), Vec::new());
    out.meta(format!("shard {}", label(shard)))?;

    let mut doms = String::new();

    let count = work.collect(shard, scratch, &mut |pair: &Pair<'_>| {
        let mut line = out.line()?;

        line.text(work.queries.name(pair.query))?;
        line.bytes(pair.target)?;
        line.bytes(pair.pass)?;

        for score in pair.scores {
            match score {
                Some(score) => line.num(*score as f64)?,
                None => line.missing()?,
            };
        }

        match pair.inc {
            Some(inc) => line.num(inc as f64)?,
            None => line.missing()?,
        };

        match pair.doms {
            [] => line.missing()?,
            scores => {
                doms.clear();
                for (i, score) in scores.iter().enumerate() {
                    let separator = match i {
                        0 => "",
                        _ => ",",
                    };
                    write!(doms, "{separator}{score:.1}").expect("a String takes what it is given");
                }

                line.text(&doms)?
            }
        };

        line.end()?;
        Ok(())
    })?;

    Ok((out.into_inner(), count))
}

/// The table's columns, which are fixed before a row is read: a stream cannot
/// widen a column once the header is out.
fn schema(meta: &Meta) -> Schema {
    let mut columns = vec![
        // a query's rows are adjacent, so they line up with each other without
        // the column being padded to the widest family name in Pfam
        Column::new("query").ragged(),
        Column::new("target").min_width(TARGET),
        Column::new("pass").min_width(meta.runs.len()),
    ];

    // named for the run rather than the tool: the
    // column header is what says which cell of the sweep a score came from
    columns.extend(
        meta.runs
            .iter()
            .map(|run| Column::new(run.name.clone()).fixed(1).min_width(SCORE)),
    );

    columns.push(Column::new("inc").min_width(3));
    columns.push(Column::new("dom").ragged());

    Schema::new(columns)
}

// -------------------------------------------------------------------- read

/// Where this table's columns sit: query, target, pass, one per run, inc, dom.
pub fn layout(meta: &Meta) -> Layout {
    let runs = meta.runs.len();

    Layout {
        fields: SCORES_AT + runs + 2,
        pass: 2,
        dom: SCORES_AT + runs + 1,
    }
}

pub struct Reader<R> {
    frame: Frame<R>,
    /// Where the first run's score column sits.
    scores: usize,
}

impl Reader<std::fs::File> {
    pub fn open(path: &Path) -> anyhow::Result<Reader<std::fs::File>> {
        Reader::wrap(Frame::open(path)?)
    }
}

impl<R: Read> Reader<R> {
    #[cfg(test)]
    pub fn new(src: R, file: &str) -> anyhow::Result<Reader<R>> {
        Reader::wrap(Frame::new(src, file)?)
    }

    fn wrap(mut frame: Frame<R>) -> anyhow::Result<Reader<R>> {
        ensure!(
            frame.format() == FORMAT,
            "{} does not open `{FORMAT}`; it is not a runs table",
            frame.file(),
        );

        frame.layout(layout(&frame.meta));

        Ok(Reader {
            frame,
            scores: SCORES_AT,
        })
    }

    pub fn meta(&self) -> &Meta {
        &self.frame.meta
    }

    /// Advance to the next row. `false` at the end of the file.
    pub fn step(&mut self) -> anyhow::Result<bool> {
        self.frame.step()
    }

    /// The frame under this reader: every column but the scores and the seed
    /// flag.
    pub fn row(&self) -> &Frame<R> {
        &self.frame
    }

    /// Whether the seeding this run replayed offered the pair.
    ///
    /// True for a run that replayed no seed list: nothing was withheld from
    /// it, which is a different answer from a seeding that looked and did not
    /// find it, and the `-` is what carries that difference.
    pub fn seeded(&self, run: usize) -> bool {
        self.frame.pass().get(run) != Some(&b'-')
    }

    /// Whether one run reported the pair at all, at any score.
    pub fn present(&self, run: usize) -> bool {
        // presence, not a cutoff: a pair that survived to be
        // scored badly was not dropped, which is the
        // distinction the stages is built on
        self.frame.filled(self.scores + run)
    }

    /// What one run scored the pair.
    #[cfg(test)]
    pub fn score(&self, run: usize) -> Option<f32> {
        super::scan::score(self.frame.field(self.scores + run))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two cells and hmmer, over a pipeline that seeded: one pair both cells
    /// found and scored differently, one seeded pair only hmmer reported, and
    /// one pair seeding never found.
    const FILE: &str = "\
#= format runs 1
#= query 3 18 120
#= target 1 3 12 40
#= seed once 4.0000
#= cutoffs /x/cutoffs.tbl c=0
#= run A2.0-B4.0 nail 3.0000 seeds=once A=2.0 B=4.0
#= run full nail 9.0000 seeds=once
#= run hmmer hmmer 2.0000
#= pass A2.0-B4.0 full hmmer
# query target           pass A2.0-B4.0 full   hmmer  inc dom
# ----- ---------------- ---- --------- ------ ------ --- ---
#= shard 1
alpha MGYP000000000001 NNH  24.0      25.0   24.0   1   24.0
beta MGYP000000000002 nnH  -         -      30.0   1   30.0
gamma MGYP000000000003 --H  -         -      28.0   1   28.0
#= end 3
";

    #[test]
    fn a_cell_keeps_its_own_score() {
        let mut reader = Reader::new(FILE.as_bytes(), "test").unwrap();

        assert!(reader.step().unwrap());

        assert_eq!(reader.score(0), Some(24.0));
        assert_eq!(reader.score(1), Some(25.0));
        assert!(reader.present(0));
    }

    #[test]
    fn the_checkpoints_are_told_apart() {
        let mut reader = Reader::new(FILE.as_bytes(), "test").unwrap();

        assert!(reader.step().unwrap());
        assert!(reader.seeded(0));
        assert!(reader.present(0));

        // seeded, and then lost between there and the table
        assert!(reader.step().unwrap());
        assert!(reader.seeded(0));
        assert!(!reader.present(0));

        // never seeded at all, which is a different loss. hmmer replayed no
        // seed list, so its own character is untouched
        assert!(reader.step().unwrap());
        assert!(!reader.seeded(0));
        assert!(!reader.present(0));
        assert!(reader.seeded(2));

        assert!(!reader.step().unwrap());
    }

    #[test]
    fn a_table_of_another_format_is_refused() {
        let other = FILE.replace(FORMAT, "#= format scores 2");
        assert!(Reader::new(other.as_bytes(), "test").is_err());
    }

    #[test]
    fn a_run_that_replayed_no_seed_list_counts_as_offered_every_pair() {
        // nothing was withheld from it, which is not the same answer as a
        // seeding that looked and did not find the pair
        let text = FILE
            .replace("#= seed once 4.0000\n", "")
            .replace(" seeds=once", "")
            .replace("--H", "nnH");

        let mut reader = Reader::new(text.as_bytes(), "test").unwrap();

        assert!(reader.step().unwrap());
        assert!(reader.seeded(0));
        assert_eq!(reader.score(1), Some(25.0));
    }
}
