use std::path::PathBuf;

use anyhow::ensure;
use bioio::split::{self, Kind};
use clap::Parser;
use pipeline::{Cmd, PipelineBuilder, Progress, Step};
use tools::{hmmsearch, mmseqs, nail};

use crate::util;

#[derive(Parser, Debug)]
pub struct Args {
    bench_dir: PathBuf,

    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    #[arg(long)]
    results: Option<PathBuf>,

    #[arg(long)]
    tmp: Option<PathBuf>,

    #[arg(long)]
    dry_run: bool,
}

const NAIL_S: [f32; 3] = [9.0, 10.0, 12.0];

const MMSEQS_S: [f32; 2] = [7.5, 12.0];
const MMSEQS_MAX_SEQS: usize = 2000;

const HMMER_CPU: usize = 2;

const EVALUE: &str = "10";

pub fn main(args: Args) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(HMMER_CPU),
        "--threads needs to be a multiple of {HMMER_CPU} (for hmmer)"
    );

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

    let query_hmm = args.bench_dir.join("queries/query.hmm");
    let query_db = args.bench_dir.join("queries/queryDB/queryDB");

    let shards = util::shards(&args.bench_dir.join("targets/"))?;

    let parts_dir = tmp_dir.join("hmmer-query");
    std::fs::remove_dir_all(&parts_dir).ok();
    let parts = split::write_splits(&query_hmm, Kind::Hmm, args.threads / HMMER_CPU, &parts_dir)?;

    let mut pl = PipelineBuilder::new().step(
        Cmd::new("mkdir")
            .flag("-p")
            .path(&results_dir)
            .path(tmp_dir.join("hmmer")),
    );

    for (idx, shard) in shards {
        let scratch = tmp_dir.join(format!("shard-{idx}"));
        let mmseqs_tdb = scratch.join("targetDB/targetDB");
        let hmmer_tmp = tmp_dir.join("hmmer");

        pl = pl
            .step(Step::serial([
                MMSEQS_S.iter().fold(
                    Cmd::new("mkdir").flag("-p").path(scratch.join("targetDB")),
                    |cmd, s| cmd.path(scratch.join(format!("alnDB-s{s:.1}"))),
                ),
                Cmd::new(&mmseqs)
                    .sub("createdb")
                    .path(&shard)
                    .path(&mmseqs_tdb),
            ]))
            .step(Step::serial(NAIL_S.iter().map(|&s| {
                Cmd::new(&nail)
                    .sub("search")
                    .arg("--mmseqs-path", &mmseqs)
                    .arg("-t", args.threads)
                    .arg("--tmp-dir", scratch.join(format!("nail-s{s:.1}")))
                    .arg("--mmseqs-s", format!("{s:.1}"))
                    .arg("--seed-mode", "prog")
                    .arg("-E", EVALUE)
                    .arg(
                        "--tbl-out",
                        results_dir.join(format!("nail-s{s:.1}.{idx}.tbl")),
                    )
                    .flag("--allow-overwrite")
                    .path(&query_hmm)
                    .path(&shard)
            })))
            .step(Step::serial(
                MMSEQS_S
                    .iter()
                    .flat_map(|&s| {
                        let adb = scratch.join(format!("alnDB-s{s:.1}/alnDB"));
                        [
                            Cmd::new(&mmseqs)
                                .sub("search")
                                .arg("--threads", args.threads)
                                .arg("-s", format!("{s:.1}"))
                                .arg("--max-seqs", MMSEQS_MAX_SEQS)
                                .arg("-e", EVALUE)
                                .path(&query_db)
                                .path(&mmseqs_tdb)
                                .path(&adb)
                                .path(scratch.join(format!("mmseqs-work-s{s:.1}"))),
                            Cmd::new(&mmseqs)
                                .sub("convertalis")
                                .arg("--format-mode", 0)
                                .path(&query_db)
                                .path(&mmseqs_tdb)
                                .path(&adb)
                                .path(
                                    results_dir.join(format!(
                                        "mmseqs-s{s:.1}-ms{MMSEQS_MAX_SEQS}.{idx}.tbl"
                                    )),
                                ),
                        ]
                    })
                    .collect::<Vec<_>>(),
            ))
            .step(Step::batch(
                parts.len(),
                parts.iter().enumerate().map(|(i, part)| {
                    Cmd::new(&hmmsearch)
                        .arg("--cpu", HMMER_CPU)
                        .arg("--tblout", hmmer_tmp.join(format!("{i}.tbl")))
                        .arg("--domtblout", hmmer_tmp.join(format!("{i}.domtbl")))
                        .arg("-E", EVALUE)
                        .path(part)
                        .path(&shard)
                }),
            ))
            .step(Step::serial([
                cat(
                    (0..parts.len()).map(|i| hmmer_tmp.join(format!("{i}.tbl"))),
                    results_dir.join(format!("hmmer.{idx}.tbl")),
                ),
                cat(
                    (0..parts.len()).map(|i| hmmer_tmp.join(format!("{i}.domtbl"))),
                    results_dir.join(format!("hmmer.{idx}.domtbl")),
                ),
            ]))
            // only this shard's scratch: the query splits and the hmmer parts
            // live on for the shards that follow
            .step(Cmd::new("rm").flag("-rf").path(&scratch));
    }

    let pipeline = pl
        .stderr_dir(tmp_dir.join("stderr"))
        .sink(Progress::new())
        .build()?;

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
