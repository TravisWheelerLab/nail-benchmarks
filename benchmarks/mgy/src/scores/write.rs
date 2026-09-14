//! Writing recall's `scores.tbl`: the shape of the table, and what a row of it
//! says.
//!
//! One column per tool rather than per run, which is recall's premise: a sweep
//! of a prefilter changes which pairs a tool reports, not what it scores them.
//! The shards themselves are collected by [`super::collect`], which serves
//! every table here and knows what none of them look like.

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, ensure};

use tabl::{Column, Schema, Stream, Widths};

use util::ledger::{self, Ledger};

use super::collect::{self, Job};
use super::shard::{Count, Pair, Scratch, Shard};
use super::{Cutoffs, Meta, Queries, Tool, label, runs, sizes, tools};

/// How wide a target name is written. MGnify's are sixteen characters, and a
/// wider one overruns rather than widening the column, since the header has
/// gone out before the first row is read.
const TARGET: usize = 16;

/// How wide a score is written: four digits, a point and a place.
const SCORE: usize = 6;

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

/// Read a finished pipeline directory into `scores.tbl`.
pub fn collect(args: Args<'_>) -> anyhow::Result<Count> {
    let ran = Ledger::load(args.dir)?;
    ledger::warn(ran.failed(), "command(s)");

    let columns = runs(&ran)?;
    ensure!(!columns.is_empty(), "no finished runs in {}", args.dir.display());

    let shards = ran.shards();
    ensure!(!shards.is_empty(), "no shards in {}", args.dir.display());

    let queries = Queries::from_hmm(args.query_hmm)?;
    let cutoffs = Cutoffs::read(args.cutoffs, args.c, &queries)?;

    let meta = Meta {
        query: queries.size,
        targets: sizes::of(args.targets, &shards)?,
        // recall never seeds, so there is no stage to time and no line
        seeds: Vec::new(),
        cutoffs: collect::absolute(args.cutoffs),
        c: args.c,
        runs: columns.iter().map(|column| column.run.clone()).collect(),
    };

    // recall holds every tool to what hmmer found, and the domain breakdown is
    // hmmer's alone, so which run is hmmer's has to be a single answer
    let hmmer = meta.hmmer()?;

    let tools = tools(&meta.runs);
    let schema = schema(meta.runs.len(), &tools);
    let widths = schema.widths();

    let file = std::fs::File::create(args.out)
        .with_context(|| format!("failed to create {}", args.out.display()))?;
    let mut out = BufWriter::with_capacity(1 << 20, file);

    meta.write(super::FORMAT, &mut out)?;
    Stream::new(schema.clone(), widths.clone(), &mut out).header()?;

    let results = args.dir.join("results");
    let work = Shard {
        results: &results,
        runs: &columns,
        queries: &queries,
        cutoffs: &cutoffs,
        hmmer: Some(hmmer),
        // recall never seeds: nail's prefilter is part of its search, so
        // there is no seed list beside its results and no stage to time
        seeds: false,
    };

    let job = Job {
        work: &work,
        shards: &shards,
        threads: args.threads,
        mem: args.mem,
    };

    let count = collect::blocks(
        job,
        |work, shard, scratch| block(work, shard, scratch, &schema, &widths, &tools),
        &mut out,
    )?;

    writeln!(out, "#= end {}", count.rows)?;
    out.flush()?;

    if count.disagreements > 0 {
        eprintln!(
            "warning: {} pair(s) were scored differently by two runs of one tool; \
             the table carries the first of them",
            count.disagreements
        );
    }

    Ok(count)
}

/// One shard's rows, rendered into memory.
fn block(
    work: &Shard<'_>,
    shard: &str,
    scratch: &mut Scratch,
    schema: &Schema,
    widths: &Widths,
    tools: &[Tool],
) -> anyhow::Result<(Vec<u8>, Count)> {
    let mut out = Stream::continued(schema.clone(), widths.clone(), Vec::new());
    out.meta(format!("shard {}", label(shard)))?;

    let mut doms = String::new();
    let mut by_tool = ByTool::default();

    let mut count = work.collect(shard, scratch, &mut |pair: &Pair<'_>| {
        by_tool.take(work.runs, pair);

        let mut line = out.line()?;

        line.text(work.queries.name(pair.query))?;
        line.bytes(pair.target)?;
        line.bytes(pair.pass)?;

        for tool in tools {
            match by_tool.scores[tool.at()] {
                Some(score) => line.num(score as f64)?,
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

    count.disagreements = by_tool.disagreements;

    Ok((out.into_inner(), count))
}

/// recall's three score columns, folded out of a pair's per-run scores.
///
/// A tool's score is whichever of its runs reported the pair first, which is
/// the run of lowest index: runs of one tool agree on a pair's score by
/// construction, so a disagreement is worth counting rather than picking
/// between. That premise is recall's -- a sweep of a prefilter changes which
/// pairs a tool reports, not what it scores them -- so the fold lives here
/// rather than in the collector, which serves tables that have no such
/// premise.
#[derive(Default)]
struct ByTool {
    scores: [Option<f32>; 3],
    disagreements: u64,
}

impl ByTool {
    fn take(&mut self, runs: &[super::Column], pair: &Pair<'_>) {
        self.scores = [None; 3];

        for (at, score) in pair.scores.iter().enumerate() {
            let Some(best) = *score else { continue };
            let slot = &mut self.scores[runs[at].run.tool.at()];

            match *slot {
                None => *slot = Some(best),
                Some(first) if first != best => self.disagreements += 1,
                Some(_) => {}
            }
        }
    }
}

/// The table's columns, which are fixed before a row is read: a stream cannot
/// widen a column once the header is out.
fn schema(runs: usize, tools: &[Tool]) -> Schema {
    let mut columns = vec![
        // a query's rows are adjacent, so they line up with each other without
        // the column being padded to the widest family name in Pfam
        Column::new("query").ragged(),
        Column::new("target").min_width(TARGET),
        Column::new("pass").min_width(runs),
    ];

    columns.extend(
        tools
            .iter()
            .map(|tool| Column::new(tool.to_string()).fixed(1).min_width(SCORE)),
    );

    columns.push(Column::new("inc").min_width(3));
    columns.push(Column::new("dom").ragged());

    Schema::new(columns)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::BTreeMap;
    use std::path::PathBuf;

    use util::ledger::Row;

    /// The runs, in the order the ledger declares them: two of nail at
    /// different sensitivities, one of mmseqs, one of hmmer.
    const RUNS: [(&str, &str); 4] = [
        ("nail-a", "nail"),
        ("nail-b", "nail"),
        ("mm", "mmseqs"),
        ("hmmer", "hmmer"),
    ];

    fn tmp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("mgy-scores-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write(path: PathBuf, text: &str) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(path, text).unwrap();
    }

    /// One row of nail's table: target query tstart tend qstart qend score
    /// bias evalue cellfrac.
    fn nail(query: &str, target: &str, score: f32) -> String {
        format!("{target} {query} 1 10 1 10 {score:.1} 0.0 1e-10 0.1\n")
    }

    /// One row of blast's twelve-column tabular format, which is what mmseqs
    /// `convertalis` writes.
    fn mmseqs(query: &str, target: &str, score: f32) -> String {
        format!("{query} {target} 31.0 100 60 2 1 10 1 10 1e-10 {score:.1}\n")
    }

    /// One row of hmmer's `--tblout`, whose eighteenth field is `inc`.
    fn hmmer(query: &str, target: &str, score: f32, inc: u32) -> String {
        format!(
            "{target} - {query} - 1e-10 {score:.1} 0.5 1e-10 {score:.1} 0.5 1.1 1 0 0 1 1 1 {inc} FL=0\n"
        )
    }

    /// One row of hmmer's `--domtblout`, one per domain.
    fn dom(query: &str, target: &str, score: f32) -> String {
        format!(
            "{target} - 100 {query} - 50 1e-10 99.9 0.5 1 1 1e-12 1e-10 {score:.1} 0.5 1 40 1 40 1 40 0.97 FL=0\n"
        )
    }

    /// A pipeline directory holding two shards of results, with a pair of
    /// every kind the table has to tell apart.
    fn pipeline(dir: &Path) {
        let results = dir.join("results");

        // shard 1: MGnify names, which key as the numbers they carry
        let (t1, t2, t3) = (
            "MGYP000000000001",
            "MGYP000000000002",
            "MGYP000000000003",
        );

        write(
            results.join("nail-a.1.tbl"),
            &format!(
                "{}{}{}{}",
                // both nail runs find this one, over alpha's cutoff
                nail("alpha", t1, 25.0),
                // under every cutoff and reported by nothing else: no row
                nail("alpha", t3, 5.0),
                // over beta's cutoff
                nail("beta", t1, 15.0),
                // reported twice, and the better of them is the score
                nail("beta", t1, 12.0),
            ),
        );

        write(results.join("nail-b.1.tbl"), &nail("alpha", t1, 25.0));

        write(
            results.join("mm.1.tbl"),
            &format!(
                "{}{}",
                mmseqs("alpha", t1, 35.0),
                // beta has no mmseqs cutoff, so this cannot pass
                mmseqs("beta", t1, 50.0),
            ),
        );

        write(
            results.join("hmmer.1.tbl"),
            &format!(
                "{}{}",
                hmmer("alpha", t1, 24.0, 1),
                // hmmer alone, and under the cutoff: a row all the same
                hmmer("alpha", t2, 5.0, 0),
            ),
        );

        write(
            results.join("hmmer.1.domtbl"),
            &format!(
                "{}{}{}",
                dom("alpha", t1, 24.0),
                dom("alpha", t2, 3.0),
                dom("alpha", t2, 2.0),
            ),
        );

        // shard 2: names that are not MGnify's, which are interned and ranked
        write(
            results.join("nail-a.2.tbl"),
            &format!("{}{}", nail("gamma", "seqA", 100.0), nail("beta", "seqB", 11.0)),
        );

        write(results.join("nail-b.2.tbl"), "");
        write(results.join("mm.2.tbl"), "");
        write(results.join("hmmer.2.tbl"), &hmmer("gamma", "seqA", 100.0, 2));
        write(
            results.join("hmmer.2.domtbl"),
            &format!("{}{}", dom("gamma", "seqA", 60.0), dom("gamma", "seqA", 40.0)),
        );

        let rows: Vec<Row> = RUNS
            .iter()
            .flat_map(|(name, tool)| {
                ["1", "2"].into_iter().map(move |shard| Row {
                    name: name.to_string(),
                    tool: tool.to_string(),
                    shard: shard.to_string(),
                    stage: String::new(),
                    params: BTreeMap::from([("s".to_string(), "9.0".to_string())]),
                    wall_s: Some(1.5),
                    cpu_s: None,
                    max_rss_kb: None,
                })
            })
            .collect();

        util::ledger::Ledger::from_rows(rows)
            .write(&dir.join("ledger.tbl"))
            .unwrap();

        write(
            dir.join("targets/1.fa"),
            ">MGYP000000000001\nAAAA\n>MGYP000000000002\nCCCC\n>MGYP000000000003\nDDDD\n",
        );
        write(dir.join("targets/2.fa"), ">seqA\nEEEE\n>seqB\nFFFF\n");

        write(
            dir.join("query.hmm"),
            &format!(
                "{}{}{}",
                util::profile("alpha", 4),
                util::profile("beta", 6),
                util::profile("gamma", 8)
            ),
        );

        write(
            dir.join("cutoffs.tbl"),
            "# family nail_1 mmseqs_1\n\
             # ------ ------ --------\n\
               alpha  20.0   30.0\n\
               beta   10.0   0.0\n\
               gamma  0.0    0.0\n",
        );
    }

    fn collected(dir: &Path, out: &Path, threads: usize) -> String {
        collect(Args {
            dir,
            query_hmm: &dir.join("query.hmm"),
            targets: &dir.join("targets"),
            cutoffs: &dir.join("cutoffs.tbl"),
            c: 0,
            out,
            threads,
            mem: 1 << 30,
        })
        .unwrap();

        std::fs::read_to_string(out).unwrap()
    }

    #[test]
    fn a_pipeline_becomes_a_table() {
        let dir = tmp("table");
        pipeline(&dir);

        let text = collected(&dir, &dir.join("scores.tbl"), 1);
        let rows: String = text
            .lines()
            .skip_while(|line| !line.starts_with("# query"))
            .map(|line| format!("{line}\n"))
            .collect();

        assert_eq!(
            rows,
            "\
# query target           pass nail   mmseqs hmmer  inc dom
# ----- ---------------- ---- ------ ------ ------ --- ---
#= shard 1
alpha MGYP000000000001 NNMH 25.0   35.0   24.0   1   24.0
alpha MGYP000000000002 nnmh -      -      5.0    0   3.0,2.0
beta MGYP000000000001 Nnmh 15.0   50.0   -      -   -
#= shard 2
beta seqB             Nnmh 11.0   -      -      -   -
gamma seqA             nnmh 100.0  -      100.0  2   60.0,40.0
#= end 5
"
        );

        // the preamble says what was searched and what it was judged by
        assert!(text.starts_with("#= format scores 2\n#= query 3 18 "));
        assert!(text.contains("\n#= target 1 3 12 "));
        assert!(text.contains("\n#= run nail-a nail 3.0000 s=9.0\n"));
        assert!(text.contains("\n#= pass nail-a nail-b mm hmmer\n"));
    }

    #[test]
    fn the_threads_do_not_show_in_the_file() {
        let dir = tmp("threads");
        pipeline(&dir);

        let one = collected(&dir, &dir.join("one.tbl"), 1);
        let three = collected(&dir, &dir.join("three.tbl"), 3);

        assert_eq!(one, three);
    }
}
