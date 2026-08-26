//! Searches every query rung against every target rung at a fixed thread count.
//!
//! The loops run ascending on both axes, so the cheap corner of the grid lands
//! first and a run that gets killed still leaves a usable surface behind.

use std::path::{Path, PathBuf};

use anyhow::{bail, ensure, Context};
use clap::Parser;

use bioio::split::{self, Kind};
use pail::{Cmd, PipelineBuilder, Progress, Step, Table};
use tools::{hmmsearch, mmseqs, nail};

// hmmsearch doesn't scale past a couple of threads, so its query gets split
// threads/HMMER_CPU ways and the parts run at the same time. without this the
// benchmark measures hmmer's thread ceiling instead of how it scales
const HMMER_CPU: usize = 2;

// every tool reports down to here, so they can be compared
const EVALUE: &str = "10";

#[derive(Parser, Debug)]
pub struct Args {
    bench_dir: PathBuf,

    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    /// How many times to time each cell. The analysis wants the minimum.
    #[arg(long, default_value_t = 1)]
    reps: usize,

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
    ensure!(args.reps > 0, "--reps needs to be at least 1");

    let nail = nail()?;
    let mmseqs = mmseqs()?;
    let hmmsearch = hmmsearch()?;

    let results_dir = match args.results {
        Some(ref dir) => dir.to_owned(),
        None => args.bench_dir.join("results"),
    };

    let tmp_dir = match args.tmp {
        Some(ref dir) => dir.to_owned(),
        None => args.bench_dir.join("tmp"),
    };

    let queries = rungs(&args.bench_dir.join("queries"), None)?;
    let targets = rungs(&args.bench_dir.join("targets"), Some("fa"))?;

    let scratch = tmp_dir.join("scratch");
    let cell = scratch.join("cell");
    let target_db = scratch.join("targetDB/targetDB");

    // split once per query rung rather than once per cell: the parts don't
    // depend on the target, and how many come back is how wide the batch is.
    // whatever a previous run left would be searched as if it belonged
    let mut parts = Vec::with_capacity(queries.len());
    for &q in &queries {
        let dir = tmp_dir.join(format!("hmmer-query/{q}"));
        std::fs::remove_dir_all(&dir).ok();
        parts.push(split::write_splits(
            args.bench_dir.join(format!("queries/{q}/query.hmm")),
            Kind::Hmm,
            args.threads / HMMER_CPU,
            &dir,
        )?);
    }

    let mut pl = PipelineBuilder::new().step(Cmd::new("mkdir").flag("-p").path(&results_dir));

    for &t in &targets {
        let target_fa = args.bench_dir.join(format!("targets/{t}.fa"));

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
                Cmd::new(&mmseqs)
                    .name("createdb")
                    .sub("createdb")
                    .arg("--threads", args.threads)
                    .path(&target_fa)
                    .path(&target_db),
            ])
            .name(format!("prep.t{t}")),
        );

        for (&q, parts) in queries.iter().zip(&parts) {
            let query_hmm = args.bench_dir.join(format!("queries/{q}/query.hmm"));
            let query_db = args.bench_dir.join(format!("queries/{q}/queryDB/queryDB"));

            for rep in 1..=args.reps {
                let at = move |cmd: Cmd, tool: &str| {
                    cmd.field("tool", tool)
                        .field("q", q)
                        .field("t", t)
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
                                .path(cell.join("mmseqs/alnDB"))
                                .path(cell.join("hmmer")),
                        ])
                        .name(format!("prep.q{q}.t{t}.r{rep}")),
                    )
                    .step(
                        Step::serial([at(
                            Cmd::new(&nail)
                                .sub("search")
                                .arg("--mmseqs-path", &mmseqs)
                                .arg("-t", args.threads)
                                .arg("--tmp-dir", cell.join("nail"))
                                .arg("-E", EVALUE)
                                .arg("--tbl-out", results_dir.join(format!("nail.q{q}.t{t}.tbl")))
                                .flag("--allow-overwrite")
                                .path(&query_hmm)
                                .path(&target_fa),
                            "nail",
                        )])
                        .name(format!("nail.q{q}.t{t}.r{rep}")),
                    )
                    .step(
                        Step::serial([at(
                            Cmd::new(&mmseqs)
                                .sub("search")
                                .arg("--threads", args.threads)
                                .arg("-e", EVALUE)
                                .path(&query_db)
                                .path(&target_db)
                                .path(cell.join("mmseqs/alnDB/alnDB"))
                                .path(cell.join("mmseqs/work")),
                            "mmseqs",
                        )])
                        .name(format!("mmseqs.q{q}.t{t}.r{rep}")),
                    )
                    .step(
                        Step::serial([Cmd::new(&mmseqs)
                            .sub("convertalis")
                            .arg("--threads", args.threads)
                            .arg("--format-mode", 0)
                            .path(&query_db)
                            .path(&target_db)
                            .path(cell.join("mmseqs/alnDB/alnDB"))
                            .path(results_dir.join(format!("mmseqs.q{q}.t{t}.tbl")))])
                        .name(format!("convert.q{q}.t{t}.r{rep}")),
                    )
                    .step(
                        Step::batched(
                            parts.len(),
                            parts.iter().enumerate().map(|(i, part)| {
                                at(
                                    Cmd::new(&hmmsearch)
                                        .name(i.to_string())
                                        .arg("--cpu", HMMER_CPU)
                                        .arg("--tblout", cell.join(format!("hmmer/{i}.tbl")))
                                        .arg("--domtblout", cell.join(format!("hmmer/{i}.domtbl")))
                                        .arg("-E", EVALUE)
                                        .path(part)
                                        .path(&target_fa),
                                    "hmmer",
                                )
                            }),
                        )
                        .name(format!("hmmer.q{q}.t{t}.r{rep}")),
                    )
                    .step(
                        Step::serial([
                            cat(
                                (0..parts.len()).map(|i| cell.join(format!("hmmer/{i}.tbl"))),
                                results_dir.join(format!("hmmer.q{q}.t{t}.tbl")),
                            )
                            .name("tbl"),
                            cat(
                                (0..parts.len()).map(|i| cell.join(format!("hmmer/{i}.domtbl"))),
                                results_dir.join(format!("hmmer.q{q}.t{t}.domtbl")),
                            )
                            .name("domtbl"),
                        ])
                        .name(format!("cat.q{q}.t{t}.r{rep}")),
                    );
            }
        }
    }

    pl = pl.step(Cmd::new("rm").name("clean").flag("-rf").path(&scratch));

    let pipeline = pl
        .stderr_dir(tmp_dir.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(args.bench_dir.join("runs.tbl")))
        .build()?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}

/// The rungs a build left behind, read off the directory rather than guessed
/// from arguments, so a run always matches the benchmark it was pointed at.
///
/// With an extension the rungs are files named `<n>.<ext>`; without one they are
/// directories named `<n>`.
fn rungs(dir: &Path, ext: Option<&str>) -> anyhow::Result<Vec<usize>> {
    let entries =
        std::fs::read_dir(dir).with_context(|| format!("failed to read {}", dir.display()))?;

    let mut out: Vec<usize> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| match ext {
            Some(ext) => p.extension().is_some_and(|x| x == ext),
            None => p.is_dir(),
        })
        .filter_map(|p| p.file_stem()?.to_str()?.parse::<usize>().ok())
        .collect();

    out.sort_unstable();

    if out.is_empty() {
        bail!("no rungs in {}", dir.display());
    }

    Ok(out)
}

/// There's no shell to expand a glob, so the parts get named one by one.
fn cat(parts: impl IntoIterator<Item = PathBuf>, into: PathBuf) -> Cmd {
    parts
        .into_iter()
        .fold(Cmd::new("cat"), |cmd, part| cmd.path(part))
        .stdout_to(into)
}
