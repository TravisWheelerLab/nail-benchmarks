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
//! #= seed <shard> <wall_s>
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

use super::collect::{self, Job};
use super::frame::{Frame, Layout};
use super::shard::{Count, Pair, Scratch, Shard};
use super::{Cutoffs, Meta, Queries, label, runs, sizes};

/// The `#= format` line a runs table opens with.
pub const FORMAT: &str = "#= format runs 1";

/// How wide a target name is written. MGnify's are sixteen characters, and a
/// wider one overruns rather than widening the column, since the header has
/// gone out before the first row is read.
const TARGET: usize = 16;

/// How wide a score is written: four digits, a point and a place.
const SCORE: usize = 6;

/// Where the score columns start: after query, target, pass and, where the
/// pipeline seeded, `seeded`.
fn scores_at(seeded: bool) -> usize {
    3 + usize::from(seeded)
}

pub struct Args<'a> {
    /// The pipeline directory: `ledger.tbl` and `results/`.
    pub dir: &'a Path,
    pub query_hmm: &'a Path,
    pub targets: &'a Path,
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

    // what the pipeline itself recorded about seeding: a stage row per shard.
    // its presence is what says the seed lists are there to be read, so a
    // pipeline that does not seed gets no column rather than an empty one
    let seeds: Vec<(String, f64)> = ran
        .stage("seed")
        .map(|row| (row.shard.clone(), row.wall_s.unwrap_or_default()))
        .collect();
    let seeded = !seeds.is_empty();

    let queries = Queries::from_hmm(args.query_hmm)?;
    let cutoffs = Cutoffs::read(args.cutoffs, args.c, &queries)?;

    let meta = Meta {
        query: queries.size,
        targets: sizes::of(args.targets, &shards)?,
        seeds,
        cutoffs: collect::absolute(args.cutoffs),
        c: args.c,
        runs: columns.iter().map(|column| column.run.clone()).collect(),
    };

    // every run here is measured against what hmmer found, and the domain
    // breakdown is hmmer's alone, so which run is hmmer's has to be a single
    // answer
    let hmmer = meta.hmmer()?;

    let schema = schema(&meta, seeded);
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
        seeds: seeded,
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

        if let Some(seeded) = pair.seeded {
            line.text(match seeded {
                true => "y",
                false => "n",
            })?;
        }

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
fn schema(meta: &Meta, seeded: bool) -> Schema {
    let mut columns = vec![
        // a query's rows are adjacent, so they line up with each other without
        // the column being padded to the widest family name in Pfam
        Column::new("query").ragged(),
        Column::new("target").min_width(TARGET),
        Column::new("pass").min_width(meta.runs.len()),
    ];

    if seeded {
        columns.push(Column::new("seeded"));
    }

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

/// Where this table's columns sit: query, target, pass, `seeded` where the
/// pipeline kept one, one per run, inc, dom.
pub fn layout(meta: &Meta) -> Layout {
    let scores = scores_at(!meta.seeds.is_empty());
    let runs = meta.runs.len();

    Layout {
        fields: scores + runs + 2,
        pass: 2,
        dom: scores + runs + 1,
    }
}

pub struct Reader<R> {
    frame: Frame<R>,
    /// Where the first run's score column sits.
    scores: usize,
    seeded: bool,
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

        // the preamble says whether the pipeline seeded, so the column is
        // placed by what the file declares rather than by counting its fields
        let seeded = !frame.meta.seeds.is_empty();
        let scores = scores_at(seeded);

        frame.layout(layout(&frame.meta));

        Ok(Reader {
            frame,
            scores,
            seeded,
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

    /// Whether seeding found this pair, `None` where the pipeline kept no
    /// seed list.
    pub fn seeded(&self) -> Option<bool> {
        match self.seeded {
            true => Some(self.frame.field(3) == b"y"),
            false => None,
        }
    }

    /// Whether one run reported the pair at all, at any score.
    pub fn present(&self, run: usize) -> bool {
        // presence, not a cutoff: a pair that survived to be
        // scored badly was not dropped, which is the
        // distinction the funnel is built on
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
#= seed 1 4.0000
#= cutoffs /x/cutoffs.tbl c=0
#= run A2.0-B4.0 nail 3.0000 A=2.0 B=4.0
#= run full nail 9.0000
#= run hmmer hmmer 2.0000
#= pass A2.0-B4.0 full hmmer
# query target           pass seeded A2.0-B4.0 full   hmmer  inc dom
# ----- ---------------- ---- ------ --------- ------ ------ --- ---
#= shard 1
alpha MGYP000000000001 NNH  y      24.0      25.0   24.0   1   24.0
beta MGYP000000000002 nnH  y      -         -      30.0   1   30.0
gamma MGYP000000000003 nnH  n      -         -      28.0   1   28.0
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
        assert_eq!(reader.seeded(), Some(true));
        assert!(reader.present(0));

        // seeded, and then lost between there and the table
        assert!(reader.step().unwrap());
        assert_eq!(reader.seeded(), Some(true));
        assert!(!reader.present(0));

        // never seeded at all, which is a different loss
        assert!(reader.step().unwrap());
        assert_eq!(reader.seeded(), Some(false));
        assert!(!reader.present(0));

        assert!(!reader.step().unwrap());
    }

    #[test]
    fn a_table_of_another_format_is_refused() {
        let other = FILE.replace(FORMAT, "#= format scores 2");
        assert!(Reader::new(other.as_bytes(), "test").is_err());
    }

    #[test]
    fn a_pipeline_that_never_seeded_has_no_column() {
        let text = FILE
            .replace("#= seed 1 4.0000\n", "")
            .replace(" pass seeded ", " pass ")
            .replace(" ---- ------ ", " ---- ")
            .replace("NNH  y      ", "NNH  ")
            .replace("nnH  y      ", "nnH  ")
            .replace("nnH  n      ", "nnH  ");

        let mut reader = Reader::new(text.as_bytes(), "test").unwrap();

        assert!(reader.step().unwrap());
        assert_eq!(reader.seeded(), None);
        assert_eq!(reader.score(1), Some(25.0));
    }
}
