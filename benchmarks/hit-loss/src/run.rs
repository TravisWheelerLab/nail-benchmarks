//! Seeds once, runs hmmer for truth, then runs nail once at its defaults.
//!
//! There is no grid here the way there is in cloud-search: this benchmark
//! isn't asking how a knob trades off, it's asking where a single default run
//! loses hits hmmer would have found. One nail invocation is the whole point.
//!
//! nail's own -E is set far above its default (10.0) rather than left alone.
//! `parse` tells "seeded but unreported" apart from "reported" by checking
//! whether a pair is in nail's .tbl at all -- and that only means what it
//! should if nothing gets cut by the final e-value gate on the way out. A
//! sky-high -E turns that gate off, so a pair missing from the .tbl can only
//! have died in cloud search or alignment, not at the door on the way out.
//!
//! -A, -B and every other pruning knob are left at nail's defaults, since
//! those are exactly the stages this benchmark is measuring the cost of.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use bioio::split::{self, Kind};
use pail::{Closure, Cmd, PipelineBuilder, Progress, Step, Table};
use tools::{hmmsearch, mmseqs, nail};

// hmmsearch doesn't scale past a couple of threads, so its query gets split
// threads/HMMER_CPU ways and the parts run at the same time
const HMMER_CPU: usize = 2;

// hmmer's own reporting threshold, for the truth set. unrelated to nail's -E,
// which this benchmark deliberately sets far above it
const HMMER_EVALUE: &str = "10";

// the seeding settings this benchmark searches with
const MMSEQS_S: &str = "12.0";
const SEED_MODE: &str = "prog";

#[derive(Parser, Debug)]
pub struct Args {
    /// A benchmark directory `hit-loss build` made.
    bench_dir: PathBuf,

    /// nail's -E, set far above its default so the final e-value gate can't
    /// be mistaken for a cloud/align filter. Only lower this to study the
    /// e-value gate itself.
    #[arg(long, default_value_t = 1e6, value_name = "X")]
    nail_evalue: f64,

    /// Threads per search, and the cores each search is pinned to.
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

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

    let nail = nail()?;
    let mmseqs = mmseqs()?;
    let hmmsearch = hmmsearch()?;

    let tmp_dir = args.tmp.unwrap_or_else(|| args.bench_dir.join("tmp"));

    let query_hmm = args.bench_dir.join("queries/query.hmm");
    let target_fa = args.bench_dir.join("targets/target.fa");
    let seeds = args.bench_dir.join("seeds");
    let truth = args.bench_dir.join("truth");
    let nail_tbl = args.bench_dir.join("nail.tbl");

    ensure!(
        query_hmm.is_file() && target_fa.is_file(),
        "{} doesn't look built; run `hit-loss build` first",
        args.bench_dir.display()
    );

    // how many ways the query splits is settled by the thread count, and
    // write_splits names the parts after their index, so the batch below can
    // be written without waiting to see what the split produced
    let jobs = args.threads / HMMER_CPU;
    let parts_dir = tmp_dir.join("hmmer-query");
    let parts: Vec<PathBuf> = (0..jobs)
        .map(|i| parts_dir.join(format!("{i}.hmm")))
        .collect();

    let pl = PipelineBuilder::new()
        .step(
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
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
                        .arg("-E", HMMER_EVALUE)
                        .path(part)
                        .path(&target_fa)
                        .field("stage", "truth")
                }),
            )
            .name("hmmer")
            // per command, not per step, so this asks for HMMER_CPU x parts,
            // which is --threads again. a machine with a smaller pool than
            // that won't fail, it will just run fewer of the parts at once
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
        )
        .step(
            Step::serial([Cmd::new(&nail)
                .sub("search")
                .arg("-t", args.threads)
                .arg("--seeds", seeds.join("seeds"))
                .arg("-E", args.nail_evalue)
                .arg("--tmp-dir", tmp_dir.join("align"))
                .arg("--tbl-out", &nail_tbl)
                .flag("--allow-overwrite")
                .path(&query_hmm)
                .path(&target_fa)
                .field("stage", "nail")])
            .name("nail")
            .cores(args.threads),
        )
        .stderr_dir(tmp_dir.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(args.bench_dir.join("runs.tbl")))
        .build()
        .context("failed to build the run")?;

    if args.dry_run {
        pl.dry_run();
        return Ok(());
    }

    pl.run()
}

/// There's no shell to expand a glob, so the parts get named one by one.
fn cat(parts: impl IntoIterator<Item = PathBuf>, into: PathBuf) -> Cmd {
    parts
        .into_iter()
        .fold(Cmd::new("cat"), |cmd, part| cmd.path(part))
        .stdout_to(into)
}
