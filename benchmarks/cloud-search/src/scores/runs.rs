//! `runs.tbl`: one row per pair, one score column per run.
//!
//! `-A` and `-B` constrain the dynamic programming, so a pair can come back
//! from a pruned cell with a lower score than the full run gave it. So the
//! score is kept per run, and a column per cell is a price this pipeline can
//! pay: it searches one unit per target source.
//!
//! ```text
//! #= format runs 3
//! #= query <count> <residues> <bytes>
//! #= target <shard> <count> <residues> <bytes>
//! #= seed <seeding> <wall_s>
//! #= cutoffs <path> c=<n>
//! #= run <name> <wall_s>[,<wall_s> ...] [k=v ...]
//! #                       A2.0
//! #                       B4.0
//! # query target           a5   ... full
//! # ----- ---------------- ---- --- ------
//! #= shard 1
//! 2-Hacid_dh_C MGYP000522683479 98.8 ... 98.8
//! #= end <rows>
//! ```
//!
//! A run's name is stacked over its column on the dashes it is built from, so
//! that a cell of the sweep is named in full without ruling fourteen columns
//! of padding across every row.
//!
//! A cell is a score, or one of two sentinels: `-` where the run's seeding
//! offered the pair and the run did not report it, and `.` where that run's
//! seeding never offered it at all. Absent is not a score: a zero or a NaN
//! would compare against a threshold and look like one, which for this table
//! is the whole question.
//!
//! A run is a cell of the sweep over every rep it ran. Its scores are rep
//! 1's, since nail scores a pair the same every time; the `#= run` line keeps
//! a wall time per rep.
//!
//! Whether a run cleared its family's cutoff is its score against nail's
//! cutoff in the file `#= cutoffs` names.

use std::io::{BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, ensure};

use toil::{Column, Schema, Stream, Widths};

use util::ledger::{self, Ledger};
use util::set::Set;

use crate::scores::collect::{self, Job};
use crate::scores::frame::{Frame, Layout, Row, Verdict};
use crate::scores::shard::{Count, Pair, Scratch, Shard};
use crate::scores::{Cutoffs, Meta, Named, Queries, runs};

/// The `#= format` line a runs table opens with.
pub const FORMAT: &str = "runs 3";

/// The format before reps, whose `#= run` lines carry one wall time and so
/// read as one rep.
const ONE_REP: &str = "runs 2";

/// How wide a target name is written. MGnify's are sixteen characters, and a
/// wider one overruns rather than widening the column, since the header has
/// gone out before the first row is read.
const TARGET: usize = 16;

/// How wide a score is written: four digits, a point and a place.
const SCORE: usize = 6;

/// Where the score columns start: after query and target.
pub const SCORES_AT: usize = 2;

/// The cell of a run whose seeding never offered the pair.
const NEVER: &str = ".";

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

    let results = args.dir.join("results");
    let columns = runs(&ran, &results)?;
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
        targets: crate::scores::target_sizes(args.set, &shards)?,
        seeds,
        cutoffs: collect::absolute(args.cutoffs),
        c: args.c,
        runs: columns.iter().map(|column| column.run.clone()).collect(),
        tools: ran.tools().to_vec(),
    };

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

    let mut head = Stream::new(schema.clone(), widths.clone(), &mut out);
    meta.write(&mut head)?;
    head.header()?;

    let work = Shard {
        results: &results,
        runs: &columns,
        queries: &queries,
        cutoffs: &cutoffs,
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

    Stream::continued(schema.clone(), widths.clone(), &mut out).meta("end", [count.rows])?;
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
    out.meta("shard", [shard])?;

    let count = work.collect(shard, scratch, &mut |pair: &Pair<'_>| {
        let mut line = out.line()?;

        line.text(work.queries.name(pair.query))?;
        line.bytes(pair.target)?;

        for (score, offered) in pair.scores.iter().zip(pair.offered) {
            match score {
                Some(score) => line.num(*score as f64)?,
                None if !offered => line.text(NEVER)?,
                None => line.missing()?,
            };
        }

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
    ];

    // named for the run: the column header is what says
    // which cell of the sweep a score came from. stacked on the dashes the name is built
    // from, so A10.0-B12.0-a1 rules three header lines
    // rather than fourteen columns of padding on every
    // row -- at 164 runs that padding was 46% of the file
    columns.extend(meta.runs.iter().map(|run| {
        Column::stacked(run.name.split('-'))
            .fixed(1)
            .min_width(SCORE)
    }));

    Schema::new(columns)
}

// -------------------------------------------------------------------- read

/// Where this table's columns sit: query, target, one per run.
pub fn layout(meta: &Meta) -> Layout {
    Layout {
        fields: SCORES_AT + meta.runs.len(),
    }
}

pub struct Reader<R: Read> {
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
            frame.format() == FORMAT || frame.format() == ONE_REP,
            "{} does not open `#= format {FORMAT}`; it is not a runs table",
            frame.file(),
        );

        let cutoffs = Named::read(&frame.meta.cutoffs, frame.meta.c)?;

        frame.layout(layout(&frame.meta));
        frame.verdict(Verdict {
            at: SCORES_AT,
            cutoffs,
        });

        Ok(Reader {
            frame,
            scores: SCORES_AT,
        })
    }

    pub fn meta(&self) -> &Meta {
        &self.frame.meta
    }

    /// Hand every row to `f` in file order, then check the trailer.
    pub fn each(
        &mut self,
        mut f: impl FnMut(&RunsRow) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let scores = self.scores;
        self.frame.each(|row| f(&RunsRow { row, scores }))
    }
}

/// One row of a runs table, read through its frame.
pub struct RunsRow<'r, 'a> {
    row: &'r Row<'a>,
    /// Where the first run's score column sits.
    scores: usize,
}

impl RunsRow<'_, '_> {
    /// The frame's row: every column but the scores and the seed flag.
    pub fn row(&self) -> &Row<'_> {
        self.row
    }

    /// Whether one run reported the pair at all, at any score.
    pub fn present(&self, run: usize) -> bool {
        super::scan::score(self.row.field(self.scores + run)).is_some()
    }

    /// What one run scored the pair.
    #[cfg(test)]
    pub fn score(&self, run: usize) -> Option<f32> {
        super::scan::score(self.row.field(self.scores + run))
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    use std::path::PathBuf;
    use std::sync::OnceLock;

    /// The cutoffs the fixture is held against, written where a reader can
    /// open them: a reader works a run's verdict out from its score rather
    /// than reading a character that already says so.
    fn cutoffs() -> &'static Path {
        // written once: the tests run together, and a second write truncates
        // the file a reader is part way through
        static PATH: OnceLock<PathBuf> = OnceLock::new();

        PATH.get_or_init(|| {
            let dir = std::env::temp_dir().join(format!("mgy-runs-{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();

            let path = dir.join("cutoffs.tbl");
            std::fs::write(
                &path,
                "# family nail_1 mmseqs_1 hmmer_1 hmmer_dom_1\n\
                 # ------ ------ -------- ------- -----------\n\
                   alpha  24.5   30.0     25.0    20.0\n\
                   beta   20.0   30.0     25.0    20.0\n\
                   zeta   0.0    30.0     25.0    20.0\n",
            )
            .unwrap();

            path
        })
    }

    /// One cell and the full run, over a pipeline that seeded: one pair both
    /// found and scored differently, one seeded pair the cell lost, one pair
    /// both cleared, and one pair in `zeta`, which has no cutoff for hmmer's
    /// best domain.
    pub fn file() -> String {
        format!(
            "\
#= format runs 3
#= query 3 18 120
#= target 1 3 12 40
#= seed once 4.0000
#= cutoffs {} c=0
#= run A2.0-B4.0 3.0000,5.0000,4.0000 seeds=once A=2.0 B=4.0
#= run full 9.0000,10.0000 seeds=once
#                       A2.0
# query target           B4.0 full
# ----- ---------------- ---- ------
#= shard 1
alpha MGYP000000000001 24.0 25.0
beta MGYP000000000002 -    30.0
beta MGYP000000000005 21.0 22.0
zeta MGYP000000000004 30.0 31.0
#= end 4
",
            cutoffs().display()
        )
    }

    /// What `f` makes of each row of `text`, in file order.
    fn rows<T>(text: &str, f: impl Fn(&RunsRow) -> T) -> Vec<T> {
        let mut reader = Reader::new(text.as_bytes(), "test").unwrap();
        let mut out = Vec::new();
        reader
            .each(|row| {
                out.push(f(row));
                Ok(())
            })
            .unwrap();
        out
    }

    #[test]
    fn a_cell_keeps_its_own_score() {
        let rows = rows(&file(), |r| (r.score(0), r.score(1), r.present(0)));
        assert_eq!(rows[0], (Some(24.0), Some(25.0), true));
    }

    #[test]
    fn a_run_clears_the_cutoff_on_its_score_alone() {
        let rows = rows(&file(), |r| (r.row().passed(0), r.row().passed(1)));

        // alpha's cutoff is 24.5, so the cell's 24.0 is under it and the
        // full run's 25.0 is over
        assert_eq!(rows[0], (false, true));

        // a cell that reported nothing cleared nothing
        assert_eq!(rows[1], (false, true));
    }

    #[test]
    fn a_table_of_another_format_is_refused() {
        let other = file().replace(&format!("#= format {FORMAT}"), "#= format runs 1");
        assert_ne!(other, file());
        assert!(Reader::new(other.as_bytes(), "test").is_err());
    }
}
