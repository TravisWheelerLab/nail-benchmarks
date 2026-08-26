//! How nail, mmseqs and hmmer scale as a search grows.
//!
//! Searches every query rung against every target rung at a fixed thread
//! count. The loops run ascending on both axes, so the cheap corner of the grid
//! lands first and a run that gets killed still leaves a usable surface behind.
//!
//! This is the one pipeline whose question is wall time rather than what was
//! found. It still writes its hit tables -- a search has to write somewhere,
//! and sending them to /dev/null would change what is being timed -- but
//! nothing parses them: a rung is its own measurement, and comparing what nail
//! found at one target rung against what hmmer found at another is not a
//! question. `manifest.tbl` and the two `sizes.tbl` files are the artifact.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use pail::{Cmd, PipelineBuilder, Progress, Step, Table};

use crate::inputs;
use crate::manifest;
use crate::search::{self, Bins, Dirs, Split};

#[derive(Parser, Debug)]
pub struct Args {
    /// Threads per search, and the cores each search is pinned to.
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    /// How many times to time each cell. The analysis wants the minimum.
    #[arg(long, default_value_t = 1)]
    reps: usize,

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
    ensure!(args.reps > 0, "--reps needs to be at least 1");

    let bins = Bins::find()?;

    let mut dirs = Dirs::new("search-size");
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let queries = inputs::ladder::query_rungs()?;
    let targets = inputs::ladder::target_rungs()?;

    // split once per query rung rather than once per cell: the parts don't
    // depend on the target, and a rung is searched against every target
    let splits: Vec<Split> = queries
        .iter()
        .map(|&q| {
            Split::new(
                inputs::ladder::query_hmm(q),
                dirs.tmp.join(format!("hmmer-query/{q}")),
                search::jobs(args.threads),
            )
        })
        .collect();

    let scratch = dirs.tmp.join("scratch");
    let cell = scratch.join("cell");
    let target_db = scratch.join("targetDB/targetDB");

    let mut pl = PipelineBuilder::new().step(dirs.mkdir());
    for split in &splits {
        pl = pl.step(split.step());
    }

    for &t in &targets {
        let shard = t.to_string();
        let target_fa = inputs::ladder::target(t);

        pl = pl.step(
            Step::serial([
                Cmd::new("rm").name("clean").flag("-rf").path(&scratch),
                Cmd::new("mkdir")
                    .name("dirs")
                    .flag("-p")
                    .path(scratch.join("targetDB")),
                // read the whole target once so the first tool to reach it
                // doesn't pay for the page cache on everyone else's behalf
                Cmd::new("cat").name("warm").path(&target_fa),
                // mmseqs takes every core it can find unless told otherwise, so
                // the setup around a search is held to the same count as the
                // search itself
                Cmd::new(&bins.mmseqs)
                    .name("createdb")
                    .sub("createdb")
                    .arg("--threads", args.threads)
                    .path(&target_fa)
                    .path(&target_db),
            ])
            .name(format!("prep.t{t}")),
        );

        for (&q, split) in queries.iter().zip(&splits) {
            let query_hmm = inputs::ladder::query_hmm(q);
            let query_db = inputs::ladder::query_db(q);

            for rep in 1..=args.reps {
                // the tool and the query rung make the column; the target rung
                // is the shard. a rep is a repeat of the same measurement, so
                // it is a field rather than a column of its own
                let named = |tool: &str| format!("{tool}.q{q}");
                let fields = |cmd: Cmd, tool: &str| {
                    cmd.field(manifest::NAME, named(tool))
                        .field(manifest::TOOL, tool)
                        .field(manifest::SHARD, &shard)
                        .field("q", q)
                        .field("rep", rep)
                };

                pl = pl
                    // cleaning at the front, so a cell that died leaves its
                    // scratch to look at and the next one still starts clean
                    .step(
                        Step::serial([
                            Cmd::new("rm").name("clean").flag("-rf").path(&cell),
                            Cmd::new("mkdir")
                                .name("dirs")
                                .flag("-p")
                                .path(cell.join("nail"))
                                .path(cell.join("mmseqs/alnDB")),
                        ])
                        .name(format!("prep.q{q}.t{t}.r{rep}")),
                    )
                    .step(
                        Step::serial([fields(
                            Cmd::new(&bins.nail)
                                .sub("search")
                                .arg("--mmseqs-path", &bins.mmseqs)
                                .arg("-t", args.threads)
                                .arg("--tmp-dir", cell.join("nail"))
                                .arg("-E", search::EVALUE)
                                .arg("--tbl-out", dirs.table(&named("nail"), &shard))
                                .flag("--allow-overwrite")
                                .path(&query_hmm)
                                .path(&target_fa),
                            "nail",
                        )])
                        .name(format!("nail.q{q}.t{t}.r{rep}"))
                        .cores(args.threads),
                    )
                    .step(
                        Step::serial(
                            search::Mmseqs {
                                bin: &bins.mmseqs,
                                query_db: &query_db,
                                target_db: &target_db,
                                aln_db: cell.join("mmseqs/alnDB/alnDB"),
                                work: cell.join("mmseqs/work"),
                                out: dirs.table(&named("mmseqs"), &shard),
                                threads: args.threads,
                                // mmseqs' own defaults: this benchmark is
                                // asking how it scales as it ships, not how a
                                // parameterization of it scales
                                s: None,
                                max_seqs: None,
                            }
                            .cmds()
                            .map(|cmd| fields(cmd, "mmseqs")),
                        )
                        .name(format!("mmseqs.q{q}.t{t}.r{rep}"))
                        .cores(args.threads),
                    );

                let hmmer = search::hmmer(
                    &bins.hmmsearch,
                    split,
                    &dirs,
                    &named("hmmer"),
                    &shard,
                    &target_fa,
                    &[("q", q.to_string()), ("rep", rep.to_string())],
                );

                // relabelled with the rep, so a repeated cell is told apart in
                // the progress output the way nail's and mmseqs' steps are
                pl = pl
                    .step(hmmer.search.name(format!("hmmer.q{q}.t{t}.r{rep}")))
                    .step(hmmer.cat.name(format!("cat.q{q}.t{t}.r{rep}")));
            }
        }
    }

    pl = pl.step(Cmd::new("rm").name("clean").flag("-rf").path(&scratch));

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(dirs.root.join("manifest.tbl")))
        .build()
        .context("failed to build the sweep")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}
