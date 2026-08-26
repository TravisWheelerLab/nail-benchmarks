//! Seeds once, then sweeps A x B off those seeds.
//!
//! Seeding and hmmer are in here rather than in a subcommand of their own
//! because they are cheap next to the grid: one seeding pass and one hmmer pass
//! against dozens of nail searches. They just have to happen first, and the
//! table marks them so the analysis can leave them out.
//!
//! -A and -B do nothing until after seeding, so the seed set is the same in
//! every cell. Paying mmseqs once and replaying it into all of them is what
//! keeps the runtime axis about cloud search rather than about the seeder.
//!
//! -a is left alone. It is a real knob and it carries most of the weight at
//! aggressive thresholds, but it is not what this benchmark is asking about.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use bioio::split::{self, Kind};
use pail::{Closure, Cmd, PipelineBuilder, Progress, Step, Table};
use tools::{hmmsearch, mmseqs, nail};

use crate::cell::{self, Cell};

// hmmsearch doesn't scale past a couple of threads, so its query gets split
// threads/HMMER_CPU ways and the parts run at the same time
const HMMER_CPU: usize = 2;

// every tool reports down to here, so they can be compared
const EVALUE: &str = "10";

// the seeding settings this benchmark searches with
const MMSEQS_S: &str = "12.0";
const SEED_MODE: &str = "prog";

#[derive(Parser, Debug)]
pub struct Args {
    /// A benchmark directory `cloud-search build` made.
    bench_dir: PathBuf,

    /// Local score pruning thresholds to sweep (nail's -A).
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "2,3,5,8,10,14,20,30,40",
        value_name = "X,X,..."
    )]
    alpha: Vec<f32>,

    /// Global score pruning thresholds to sweep (nail's -B).
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "4,6,8,12,16,22,32,48,64",
        value_name = "X,X,..."
    )]
    beta: Vec<f32>,

    /// Threads per search, and the cores each search is pinned to.
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    #[arg(long)]
    results: Option<PathBuf>,

    #[arg(long)]
    tmp: Option<PathBuf>,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(HMMER_CPU),
        "--threads needs to be a multiple of {HMMER_CPU} (for hmmer)"
    );
    ensure!(!args.alpha.is_empty(), "--alpha needs at least one value");
    ensure!(!args.beta.is_empty(), "--beta needs at least one value");

    let nail = nail()?;
    let mmseqs = mmseqs()?;
    let hmmsearch = hmmsearch()?;

    let results_dir = args
        .results
        .unwrap_or_else(|| args.bench_dir.join("results"));
    let tmp_dir = args.tmp.unwrap_or_else(|| args.bench_dir.join("tmp"));

    let query_hmm = args.bench_dir.join("queries/query.hmm");
    let target_fa = args.bench_dir.join("targets/target.fa");
    let seeds = args.bench_dir.join("seeds");
    let truth = args.bench_dir.join("truth");

    ensure!(
        query_hmm.is_file() && target_fa.is_file(),
        "{} doesn't look built; run `cloud-search build` first",
        args.bench_dir.display()
    );

    // how many ways the query splits is settled by the thread count, and
    // write_splits names the parts after their index, so the batch below can be
    // written without waiting to see what the split produced
    let jobs = args.threads / HMMER_CPU;
    let parts_dir = tmp_dir.join("hmmer-query");
    let parts: Vec<PathBuf> = (0..jobs)
        .map(|i| parts_dir.join(format!("{i}.hmm")))
        .collect();

    let cells = cell::cells(&args.alpha, &args.beta);

    let pl = PipelineBuilder::new()
        .step(
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(&results_dir)
                .path(&seeds)
                .path(&truth)
                // hmmsearch won't make this itself; it just fails to open its
                // output
                .path(tmp_dir.join("hmmer")),
        )
        // rust in place of a command, so it is a closure step. whatever a
        // previous run left in there would be searched as if it belonged
        .step(
            Step::from_closures([Closure::new("split", {
                let (query_hmm, parts_dir) = (query_hmm.clone(), parts_dir.clone());

                move || {
                    std::fs::remove_dir_all(&parts_dir).ok();
                    let written = split::write_splits(&query_hmm, Kind::Hmm, jobs, &parts_dir)?;

                    // empty bins are skipped, so a query with fewer models than
                    // parts comes back short and the batch below would be
                    // pointed at files that were never written
                    anyhow::ensure!(
                        written.len() == jobs,
                        "split {} into {} parts, expected {jobs}",
                        query_hmm.display(),
                        written.len()
                    );

                    Ok(())
                }
            })])
            .name("split"),
        )
        .step(
            Step::serial([Cmd::new(&nail)
                .sub("search")
                .arg("--mmseqs-path", &mmseqs)
                .arg("-t", args.threads)
                .arg("--tmp-dir", tmp_dir.join("seeding"))
                .arg("--mmseqs-s", MMSEQS_S)
                .arg("--seed-mode", SEED_MODE)
                .arg("--seeds-out", seeds.join("seeds"))
                .flag("--only-seed")
                .flag("--allow-overwrite")
                .path(&query_hmm)
                .path(&target_fa)
                .field("stage", "seed")])
            .name("seeds")
            .cores(args.threads),
        )
        .step(
            Step::batched(
                parts.len(),
                parts.iter().enumerate().map(|(i, part)| {
                    Cmd::new(&hmmsearch)
                        .name(i.to_string())
                        .arg("--cpu", HMMER_CPU)
                        .arg("--tblout", tmp_dir.join(format!("hmmer/{i}.tbl")))
                        .arg("--domtblout", tmp_dir.join(format!("hmmer/{i}.domtbl")))
                        .arg("-E", EVALUE)
                        .path(part)
                        .path(&target_fa)
                        .field("stage", "truth")
                }),
            )
            .name("hmmer")
            // per command, not per step, so this asks for HMMER_CPU x parts,
            // which is --threads again. a machine with a smaller pool than that
            // won't fail, it will just run fewer of the parts at once
            .cores(HMMER_CPU),
        )
        .step(
            Step::serial([
                cat(
                    (0..parts.len()).map(|i| tmp_dir.join(format!("hmmer/{i}.tbl"))),
                    truth.join("hmmer.tbl"),
                )
                .name("tbl"),
                cat(
                    (0..parts.len()).map(|i| tmp_dir.join(format!("hmmer/{i}.domtbl"))),
                    truth.join("hmmer.domtbl"),
                )
                .name("domtbl"),
            ])
            .name("truth"),
        );

    // the cells run in the order the grid gives them, which is ascending, so
    // the cheap corner lands first and a sweep that gets killed still leaves a
    // usable surface behind
    let pl = cells.iter().fold(pl, |pl, cell| {
        let label = cell.label();

        let cmd = Cmd::new(&nail)
            .sub("search")
            .arg("-t", args.threads)
            .arg("--seeds", seeds.join("seeds"))
            .arg("-E", EVALUE)
            .arg("--tmp-dir", tmp_dir.join("cell"))
            .arg("--tbl-out", results_dir.join(format!("{label}.tbl")))
            .flag("--allow-overwrite")
            .field("stage", "cell");

        let cmd = match cell {
            Cell::Pruned { a, b } => cmd
                .arg("-A", *a)
                .arg("-B", *b)
                .field("A", *a)
                .field("B", *b),
            Cell::Full => cmd.flag("--full-dp"),
        };

        pl.step(
            Step::serial([cmd.path(&query_hmm).path(&target_fa)])
                .name(&label)
                // every cell on the same cores, so the only thing moving
                // between them is -A and -B. without this a cell is timed
                // against whatever the scheduler felt like that second, and
                // the differences here are small enough for that to show
                .cores(args.threads),
        )
    });

    let pipeline = pl
        .stderr_dir(tmp_dir.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(args.bench_dir.join("runs.tbl")))
        .build()
        .context("failed to build the sweep")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}

/// There's no shell to expand a glob, so the parts get named one by one.
fn cat(parts: impl IntoIterator<Item = PathBuf>, into: PathBuf) -> Cmd {
    parts
        .into_iter()
        .fold(Cmd::new("cat"), |cmd, part| cmd.path(part))
        .stdout_to(into)
}
