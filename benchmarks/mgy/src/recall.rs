//! How much of what hmmer finds nail and mmseqs find, as their prefilter
//! sensitivity moves.
//!
//! This is the one pipeline that searches the whole target set rather than a
//! single shard, and the one that runs mmseqs. A shard is a unit of work
//! rather than a variable: every tool sees all of them, and a run's column in
//! the scores table spans the lot.
//!
//! hmmer is a column like nail and mmseqs. What makes it the thing the others
//! are measured against is the analysis holding them to it, not anything about
//! how it is run or recorded.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use pail::{Cmd, PipelineBuilder, Progress, Step, Table};

use crate::inputs;
use crate::search::{self, Bins, Dirs, Split};
use bench::manifest;

/// mmseqs' own default is 300, which loses hits nail's seeding keeps. 2000 is
/// what the comparison has always been run at.
const MMSEQS_MAX_SEQS: usize = 2000;

/// nail's seeding mode, matching what cloud-search and hit-loss seed with.
const SEED_MODE: &str = "prog";

/// The column hmmer's run becomes, which the other two are measured against.
const HMMER: &str = "hmmer";

#[derive(Parser, Debug)]
pub struct Args {
    /// Search only the first N shards. Every shard by default.
    #[arg(short = 'n', long, value_name = "N")]
    shards: Option<usize>,

    /// nail's --mmseqs-s values to sweep.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "9.0,10.0,12.0",
        value_name = "X,X,..."
    )]
    nail_s: Vec<f32>,

    /// mmseqs' -s values to sweep.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "7.5,12.0",
        value_name = "X,X,..."
    )]
    mmseqs_s: Vec<f32>,

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
    ensure!(!args.nail_s.is_empty(), "--nail-s needs at least one value");
    ensure!(
        !args.mmseqs_s.is_empty(),
        "--mmseqs-s needs at least one value"
    );

    let bins = Bins::find()?;

    let mut dirs = Dirs::new("recall");
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let query_hmm = inputs::fixed::query_hmm();
    let query_db = inputs::fixed::query_db();

    ensure!(
        query_hmm.is_file(),
        "{} is missing; run `mgy build` first",
        query_hmm.display()
    );

    let mut shards = inputs::fixed::shards()?;
    if let Some(n) = args.shards {
        ensure!(n > 0, "-n 0 would leave nothing to search");
        if n > shards.len() {
            eprintln!(
                "warning: asked for {n} shards but there are only {}; using all of them",
                shards.len()
            );
        }
        shards.truncate(n);
    }

    // split once rather than once per shard: the parts don't depend on the
    // target, and whatever a previous run left would be searched as if it
    // belonged
    let split = Split::new(
        &query_hmm,
        dirs.tmp.join("hmmer-query"),
        search::jobs(args.threads),
    );

    let mut pl = PipelineBuilder::new().step(dirs.mkdir()).step(split.step());

    for (idx, target) in &shards {
        let shard = idx.to_string();
        let scratch = dirs.tmp.join(format!("shard-{idx}"));
        let target_db = scratch.join("targetDB/targetDB");

        pl = pl.step(
            Step::serial([
                args.mmseqs_s.iter().fold(
                    Cmd::new("mkdir")
                        .name("dirs")
                        .flag("-p")
                        .path(scratch.join("targetDB")),
                    |cmd, s| cmd.path(scratch.join(format!("alnDB-s{s:.1}"))),
                ),
                // mmseqs takes every core it can find unless told otherwise, so
                // the setup around a search is held to the same count as the
                // search itself
                Cmd::new(&bins.mmseqs)
                    .name("createdb")
                    .sub("createdb")
                    .arg("--threads", args.threads)
                    .path(target)
                    .path(&target_db),
            ])
            .name(format!("prep.{idx}")),
        );

        // ---- nail

        pl = pl.step(
            Step::serial(args.nail_s.iter().map(|&s| {
                let name = format!("nail-s{s:.1}");

                Cmd::new(&bins.nail)
                    .sub("search")
                    .arg("--mmseqs-path", &bins.mmseqs)
                    .arg("-t", args.threads)
                    .arg("--tmp-dir", scratch.join(&name))
                    .arg("--mmseqs-s", format!("{s:.1}"))
                    .arg("--seed-mode", SEED_MODE)
                    .arg("-E", search::EVALUE)
                    .arg("--tbl-out", dirs.table(&name, &shard))
                    .flag("--allow-overwrite")
                    .path(&query_hmm)
                    .path(target)
                    .field(manifest::NAME, &name)
                    .field(manifest::TOOL, "nail")
                    .field(manifest::SHARD, &shard)
                    .field("s", format!("{s:.1}"))
            }))
            .name(format!("nail.{idx}"))
            .cores(args.threads),
        );

        // ---- mmseqs
        //
        // the search does the work and the conversion writes the table, so both
        // carry the run's fields: what the column cost is the two together.

        pl = pl.step(
            Step::serial(
                args.mmseqs_s
                    .iter()
                    .flat_map(|&s| {
                        let name = format!("mmseqs-s{s:.1}-ms{MMSEQS_MAX_SEQS}");

                        let cmds = search::Mmseqs {
                            bin: &bins.mmseqs,
                            query_db: &query_db,
                            target_db: &target_db,
                            aln_db: scratch.join(format!("alnDB-s{s:.1}/alnDB")),
                            work: scratch.join(format!("work-s{s:.1}")),
                            out: dirs.table(&name, &shard),
                            threads: args.threads,
                            s: Some(format!("{s:.1}")),
                            max_seqs: Some(MMSEQS_MAX_SEQS),
                        }
                        .cmds();

                        cmds.map(|cmd| {
                            cmd.field(manifest::NAME, &name)
                                .field(manifest::TOOL, "mmseqs")
                                .field(manifest::SHARD, &shard)
                                .field("s", format!("{s:.1}"))
                        })
                    })
                    .collect::<Vec<_>>(),
            )
            .name(format!("mmseqs.{idx}"))
            .cores(args.threads),
        );

        // ---- hmmer

        let hmmer = search::hmmer(&bins.hmmsearch, &split, &dirs, HMMER, &shard, target, &[]);
        pl = pl.step(hmmer.search).step(hmmer.cat);

        // only this shard's scratch: the query splits live on for the shards
        // that follow
        pl = pl.step(Cmd::new("rm").name("clean").flag("-rf").path(&scratch));
    }

    let pipeline = pl
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
