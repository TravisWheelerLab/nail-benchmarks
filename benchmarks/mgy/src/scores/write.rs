//! Writing recall's `scores.tbl`: the shape of the table, and the shards that
//! fill it.
//!
//! Every target lives in exactly one shard, so a shard's rows are a run of the
//! file that nothing outside it belongs in. Shards are collected in parallel
//! and written in order, which is what gets a sort of four billion rows for
//! the price of sorting each shard's own.
//!
//! The same bytes come out however many threads there were: a worker renders
//! its block into memory and the file takes them in shard order.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex, mpsc};

use anyhow::{Context, ensure};

use tabl::{Column, Schema, Stream, Widths};

use util::ledger::{self, Ledger};
use util::manifest;

use super::shard::{Count, Pair, Scratch, Shard};
use super::{Cutoffs, Meta, Queries, Tool, label, runs, sizes, tools};

/// How wide a target name is written. MGnify's are sixteen characters, and a
/// wider one overruns rather than widening the column, since the header has
/// gone out before the first row is read.
const TARGET: usize = 16;

/// How wide a score is written: four digits, a point and a place.
const SCORE: usize = 6;

/// What a shard costs to hold while it is collected, beside the share of its
/// input that [`estimate`] allows for.
const OVERHEAD: u64 = 64 << 20;

/// How much of a shard's input it holds at the peak, in eighths.
///
/// Measured rather than reasoned: a synthetic shard of 4.0 GB peaked at 3.1 GB
/// resident with one worker, which is the hits and the domains it was read
/// into, the block it was rendered to, and what the allocator was holding
/// while each of those doubled.
const PEAK: u64 = 6;

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
        cutoffs: absolute(args.cutoffs),
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

    meta.write(&mut out)?;
    Stream::new(schema.clone(), widths.clone(), &mut out).header()?;

    let results = args.dir.join("results");
    let work = Shard {
        results: &results,
        runs: &columns,
        queries: &queries,
        cutoffs: &cutoffs,
        hmmer: Some(hmmer),
    };

    let count = blocks(&work, &shards, &args, &schema, &widths, &tools, &mut out)?;

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

/// Collect every shard and write the blocks in order.
fn blocks(
    work: &Shard<'_>,
    shards: &[String],
    args: &Args<'_>,
    schema: &Schema,
    widths: &Widths,
    tools: &[Tool],
    out: &mut impl Write,
) -> anyhow::Result<Count> {
    let ticket = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    let budget = Budget::new(args.mem);
    let (send, recv) = mpsc::channel::<Message>();

    let threads = args.threads.max(1).min(shards.len());

    std::thread::scope(|scope| -> anyhow::Result<Count> {
        for _ in 0..threads {
            let send = send.clone();
            let (ticket, stop, budget) = (&ticket, &stop, &budget);

            scope.spawn(move || {
                let mut scratch = Scratch::default();

                loop {
                    if stop.load(Ordering::Relaxed) {
                        return;
                    }

                    let at = ticket.fetch_add(1, Ordering::Relaxed);
                    let Some(shard) = shards.get(at) else {
                        return;
                    };

                    let held = estimate(work, shard);
                    budget.take(at, held);

                    let message = match block(work, shard, &mut scratch, schema, widths, tools) {
                        Ok((block, count)) => Message::Block {
                            at,
                            block,
                            count,
                            held,
                        },
                        Err(error) => {
                            budget.give(held);
                            stop.store(true, Ordering::Relaxed);
                            Message::Failed {
                                shard: shard.clone(),
                                error,
                            }
                        }
                    };

                    if send.send(message).is_err() {
                        return;
                    }
                }
            });
        }

        // the senders the workers hold are clones; this one would keep the
        // channel open after they are all done
        drop(send);

        let mut pending: BTreeMap<usize, (Vec<u8>, u64)> = BTreeMap::new();
        let mut next = 0usize;
        let mut total = Count::default();
        let mut failure: Option<(String, anyhow::Error)> = None;

        for message in recv {
            match message {
                Message::Failed { shard, error } => {
                    failure.get_or_insert((shard, error));

                    // nothing more will be written, so the blocks waiting for
                    // a shard that failed give their room back rather than
                    // holding the workers still reading
                    for (_, (_, held)) in std::mem::take(&mut pending) {
                        budget.give(held);
                    }
                }
                Message::Block {
                    at,
                    block,
                    count,
                    held,
                } => {
                    if failure.is_some() {
                        budget.give(held);
                        continue;
                    }

                    total.add(count);
                    pending.insert(at, (block, held));

                    while let Some((block, held)) = pending.remove(&next) {
                        out.write_all(&block)?;
                        next += 1;

                        budget.waiting_for(next);
                        budget.give(held);
                    }
                }
            }
        }

        if let Some((shard, error)) = failure {
            return Err(error.context(format!("shard {shard}")));
        }

        Ok(total)
    })
}

enum Message {
    Block {
        at: usize,
        block: Vec<u8>,
        count: Count,
        held: u64,
    },
    Failed {
        shard: String,
        error: anyhow::Error,
    },
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

    let count = work.collect(shard, scratch, &mut |pair: &Pair<'_>| {
        let mut line = out.line()?;

        line.text(work.queries.name(pair.query))?;
        line.bytes(pair.target)?;
        line.bytes(pair.pass)?;

        for tool in tools {
            match pair.scores[tool.at()] {
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

    Ok((out.into_inner(), count))
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

/// What a shard will take to collect: its tables, the hits they become, and
/// the block they are rendered into.
fn estimate(work: &Shard<'_>, shard: &str) -> u64 {
    let bytes: u64 = work
        .runs
        .iter()
        .filter(|column| column.shards.iter().any(|covered| covered == shard))
        .map(|column| {
            let table = manifest::table_path(work.results, &column.run.name, shard);
            let dom = manifest::dom_path(work.results, &column.run.name, shard);

            [table, dom]
                .iter()
                .filter_map(|path| std::fs::metadata(path).ok())
                .map(|meta| meta.len())
                .sum::<u64>()
        })
        .sum();

    bytes * PEAK / 8 + OVERHEAD
}

/// How much memory the collectors may hold between them.
///
/// A shard is admitted when what it asks for fits under the limit, and holds
/// its share until its block has been written. One that asks for more than the
/// whole limit runs alone rather than not at all.
///
/// The shard the file is waiting for is admitted whatever is held. Without
/// that the budget can fill with blocks that have been collected and cannot be
/// written -- the file wants shard 4, the room is held by the blocks for 5, 6
/// and 7, and the worker that claimed 4 is waiting for room that only writing
/// 4 would free.
struct Budget {
    limit: u64,
    state: Mutex<State>,
    room: Condvar,
}

struct State {
    used: u64,
    /// Which shard the file is waiting for.
    next: usize,
}

impl Budget {
    fn new(limit: u64) -> Budget {
        Budget {
            limit: limit.max(1),
            state: Mutex::new(State { used: 0, next: 0 }),
            room: Condvar::new(),
        }
    }

    fn take(&self, at: usize, want: u64) {
        let mut state = self.state.lock().expect("the budget outlives its holders");

        while state.used > 0 && state.used + want > self.limit && at != state.next {
            state = self
                .room
                .wait(state)
                .expect("the budget outlives its holders");
        }

        state.used += want;
    }

    fn give(&self, want: u64) {
        let mut state = self.state.lock().expect("the budget outlives its holders");
        state.used = state.used.saturating_sub(want);

        self.room.notify_all();
    }

    /// Say which shard the file is waiting for, which is the one that cannot
    /// be made to wait.
    fn waiting_for(&self, at: usize) {
        let mut state = self.state.lock().expect("the budget outlives its holders");
        state.next = at;

        self.room.notify_all();
    }
}

/// Half of what the machine has, which is what a collector is allowed by
/// default. Falls back to eight gigabytes where the machine will not say.
pub fn ram() -> u64 {
    const FALLBACK: u64 = 8 << 30;

    #[cfg(target_os = "linux")]
    let total = std::fs::read_to_string("/proc/meminfo").ok().and_then(|text| {
        text.lines()
            .find_map(|line| line.strip_prefix("MemTotal:"))
            .and_then(|rest| rest.split_whitespace().next())
            .and_then(|kb| kb.parse::<u64>().ok())
            .map(|kb| kb * 1024)
    });

    #[cfg(target_os = "macos")]
    let total = std::process::Command::new("sysctl")
        .args(["-n", "hw.memsize"])
        .output()
        .ok()
        .and_then(|out| String::from_utf8(out.stdout).ok())
        .and_then(|text| text.trim().parse::<u64>().ok());

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    let total: Option<u64> = None;

    total.unwrap_or(FALLBACK)
}

/// A path as the `#= cutoffs` line should carry it: what was given, made
/// absolute, so a table read elsewhere still names the file it was judged by.
pub fn absolute(path: &Path) -> PathBuf {
    match path.is_absolute() {
        true => path.to_path_buf(),
        false => std::env::current_dir()
            .map(|dir| dir.join(path))
            .unwrap_or_else(|_| path.to_path_buf()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::BTreeMap;

    use util::ledger::Row;

    #[test]
    fn a_budget_admits_one_shard_bigger_than_all_of_it() {
        let budget = Budget::new(1 << 20);

        // nothing is held, so an oversized shard goes ahead rather than
        // waiting for room that will never come
        budget.take(0, 1 << 30);
        budget.give(1 << 30);
    }

    /// The shard the file is waiting for goes ahead even with the room full,
    /// since everything held is waiting on it.
    #[test]
    fn a_budget_never_holds_up_the_shard_being_written() {
        let budget = std::sync::Arc::new(Budget::new(1 << 20));
        budget.take(5, 1 << 20);
        budget.waiting_for(4);

        let waiting = budget.clone();
        let admitted = std::thread::spawn(move || waiting.take(4, 1 << 20));

        // it has to come back without anything being given back first
        admitted.join().unwrap();
    }

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
