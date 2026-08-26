//! Which stage of the nail pipeline drops the hits hmmer finds.
//!
//! Seeds once, runs hmmer, then runs nail once at its defaults.
//! There is no grid here the way there is in cloud-search: this isn't asking
//! how a knob trades off, it's asking where a single default run loses hits
//! hmmer would have found. One nail invocation is the whole point.
//!
//! nail's own -E is set far above its default (10.0) rather than left alone.
//! `parse` tells "seeded but unreported" apart from "reported" by whether a
//! pair is in nail's table at all -- and that only means what it should if
//! nothing gets cut by the final e-value gate on the way out. A sky-high -E
//! turns that gate off, so a pair missing from the table can only have died in
//! cloud search or alignment, not at the door on the way out.
//!
//! -A, -B and every other pruning knob are left at nail's defaults, since
//! those are exactly the stages this benchmark is measuring the cost of.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use pail::{Cmd, PipelineBuilder, Progress, Step, Table};

use crate::inputs;
use crate::manifest;
use crate::search::{self, Bins, Dirs, Split};

/// The seeding settings this benchmark searches with, matching cloud-search's
/// so the two are asking about the same seed set.
const MMSEQS_S: &str = "12.0";
const SEED_MODE: &str = "prog";

/// The columns the two runs become.
const RUN: &str = "nail";
const HMMER: &str = "hmmer";

#[derive(Parser, Debug)]
pub struct Args {
    /// Which target shard to search.
    #[arg(long, default_value = "1", value_name = "N")]
    shard: String,

    /// nail's -E, set far above its default so the final e-value gate can't be
    /// mistaken for a cloud/align filter. Only lower this to study the e-value
    /// gate itself.
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
        args.threads.is_multiple_of(search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        search::HMMER_CPU
    );

    let bins = Bins::find()?;

    let mut dirs = Dirs::new("hit-loss");
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let query_hmm = inputs::fixed::query_hmm();
    let target = inputs::fixed::shard(&args.shard);

    ensure!(
        query_hmm.is_file() && target.is_file(),
        "{} or {} is missing; run `mgy build` first",
        query_hmm.display(),
        target.display()
    );

    let split = Split::new(
        &query_hmm,
        dirs.tmp.join("hmmer-query"),
        search::jobs(args.threads),
    );

    let mut pl = PipelineBuilder::new()
        .step(dirs.mkdir())
        .step(split.step())
        .step(search::seed(
            &bins.nail,
            &bins.mmseqs,
            &query_hmm,
            &target,
            &args.shard,
            &dirs,
            args.threads,
            MMSEQS_S,
            SEED_MODE,
        ));

    let hmmer = search::hmmer(
        &bins.hmmsearch,
        &split,
        &dirs,
        HMMER,
        &args.shard,
        &target,
        &[],
    );
    pl = pl.step(hmmer.search).step(hmmer.cat);

    let pipeline = pl
        .step(
            Step::serial([Cmd::new(&bins.nail)
                .sub("search")
                .arg("-t", args.threads)
                .arg("--seeds", dirs.seeds(&args.shard))
                .arg("-E", args.nail_evalue)
                .arg("--tmp-dir", dirs.tmp.join("align"))
                .arg("--tbl-out", dirs.table(RUN, &args.shard))
                .flag("--allow-overwrite")
                .path(&query_hmm)
                .path(&target)
                .field(manifest::NAME, RUN)
                .field(manifest::TOOL, "nail")
                .field(manifest::SHARD, &args.shard)
                .field("E", args.nail_evalue)])
            .name(RUN)
            .cores(args.threads),
        )
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(dirs.root.join("manifest.tbl")))
        .build()
        .context("failed to build the run")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}
