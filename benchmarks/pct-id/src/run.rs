//! Runs every tool against the one target/query pair `pct-id build` made.
//!
//! nail and mmseqs sweep their prefilter sensitivity (`-s`); every other knob
//! is fixed. Every tool but last and diamond runs in both profile mode
//! (against `query.hmm`/`query.sto`) and sequence mode (against `query.fa`),
//! since the whole point of this benchmark is comparing the two.
//!
//! blast is deliberately left at its default E-value: raising it to match the
//! others makes it dramatically slower for no extra recall.

use std::path::PathBuf;

use anyhow::{bail, ensure, Context};
use clap::Parser;

use pail::{Cmd, PipelineBuilder, Progress, Step, Table};
use tools::{
    blastp, diamond, hmmsearch, lastal, lastdb, makeblastdb, mmseqs, nail, phmmer, psiblast,
};

// hmmsearch/phmmer scale poorly past a couple of threads, so the query set is
// split threads/HMMER_CPU ways and the parts run at the same time
const HMMER_CPU: usize = 2;

// every tool but blast reports down to here, so they are comparable
const EVALUE: &str = "1e9";
const MAX_SEQS: &str = "2000";

#[derive(Parser, Debug)]
pub struct Args {
    /// Which benchmark to run, naming benchmark-<size>/.
    #[arg(short, long, default_value = "toy")]
    pub size: String,

    /// nail's --mmseqs-s values to sweep.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "5.7,7.5,10.0,12.0,14.0",
        value_name = "X,X,..."
    )]
    pub nail_s: Vec<f32>,

    /// mmseqs' -s values to sweep.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "5.7,7.5,10.0,12.0,14.0",
        value_name = "X,X,..."
    )]
    pub mmseqs_s: Vec<f32>,

    #[arg(short, long, default_value_t = 24)]
    pub threads: usize,

    #[arg(long)]
    pub results: Option<PathBuf>,

    #[arg(long)]
    pub tmp: Option<PathBuf>,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(HMMER_CPU),
        "--threads needs to be a multiple of {HMMER_CPU} (for hmmer)"
    );
    ensure!(!args.nail_s.is_empty(), "--nail-s needs at least one value");
    ensure!(
        !args.mmseqs_s.is_empty(),
        "--mmseqs-s needs at least one value"
    );

    let dir = crate::dir();
    let bench_dir = dir.join(format!("benchmark-{}", args.size));
    if !bench_dir.is_dir() {
        bail!(
            "benchmark directory {} does not exist; run `pct-id build --size {}` first",
            bench_dir.display(),
            args.size
        );
    }
    // absolute, so the argv column of runs.tbl can be pasted into a shell
    let bench_dir = bench_dir.canonicalize()?;

    let results = args
        .results
        .unwrap_or_else(|| bench_dir.join("results"));
    let tmp = args.tmp.unwrap_or_else(|| bench_dir.join("tmp"));

    // a rerun's stale results and scratch would otherwise sit alongside the
    // new ones, and mmseqs refuses to overwrite an existing alignment db
    if !args.dry_run {
        std::fs::remove_dir_all(&results).ok();
        std::fs::remove_dir_all(&tmp).ok();
    }

    let query_hmm = bench_dir.join("query.hmm");
    let query_sto = bench_dir.join("query.sto");
    let query_fa = bench_dir.join("query.fa");
    let target_fa = bench_dir.join("target.fa");
    let afa_dir = bench_dir.join("afa");

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;
    let hmmsearch_bin = hmmsearch()?;
    let phmmer_bin = phmmer()?;
    let blastp_bin = blastp()?;
    let psiblast_bin = psiblast()?;
    let makeblastdb_bin = makeblastdb()?;
    let lastal_bin = lastal()?;
    let lastdb_bin = lastdb()?;
    let diamond_bin = diamond()?;

    let mmseqs_dir = tmp.join("mmseqs");
    let target_db = mmseqs_dir.join("targetDB/targetDB");
    let query_prf_db = mmseqs_dir.join("queryDB-prf/queryDB");
    let query_seq_db = mmseqs_dir.join("queryDB-seq/queryDB");
    let msa_db = mmseqs_dir.join("msaDB/msaDB");
    let blast_db = tmp.join("blast/db");
    let last_db = tmp.join("last/db");
    let diamond_db = tmp.join("diamond/db");

    // how many ways the query splits is settled by the thread count, and
    // write_splits names the parts after their index, so the batch below can
    // be written without waiting to see what the split produced
    let hmmer_jobs = args.threads / HMMER_CPU;

    let mut afa: Vec<PathBuf> = std::fs::read_dir(&afa_dir)
        .with_context(|| format!("failed to read {}", afa_dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "afa"))
        .collect();
    afa.sort();
    if afa.is_empty() {
        bail!("no .afa files in {}", afa_dir.display());
    }

    let mut pl = PipelineBuilder::new()
        .step(
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(&results)
                .path(target_db.parent().unwrap())
                .path(query_prf_db.parent().unwrap())
                .path(query_seq_db.parent().unwrap())
                .path(msa_db.parent().unwrap())
                .path(blast_db.parent().unwrap())
                .path(last_db.parent().unwrap())
                .path(diamond_db.parent().unwrap()),
        )
        .step(
            Step::serial([
                Cmd::new(&mmseqs_bin)
                    .name("createdb-target")
                    .sub("createdb")
                    .path(&target_fa)
                    .path(&target_db),
                Cmd::new(&mmseqs_bin)
                    .name("createdb-query-seq")
                    .sub("createdb")
                    .path(&query_fa)
                    .path(&query_seq_db),
                Cmd::new(&mmseqs_bin)
                    .name("convertmsa")
                    .sub("convertmsa")
                    .path(&query_sto)
                    .path(&msa_db)
                    .arg("--identifier-field", 0),
                Cmd::new(&mmseqs_bin)
                    .name("msa2profile")
                    .sub("msa2profile")
                    .path(&msa_db)
                    .path(&query_prf_db)
                    .arg("--match-mode", 1),
                Cmd::new(&makeblastdb_bin)
                    .name("makeblastdb")
                    .arg("-in", &target_fa)
                    .arg("-dbtype", "prot")
                    .arg("-out", &blast_db),
                Cmd::new(&lastdb_bin)
                    .name("lastdb")
                    .arg("-p", &last_db)
                    .path(&target_fa),
                Cmd::new(&diamond_bin)
                    .name("diamond-makedb")
                    .sub("makedb")
                    .arg("--in", &target_fa)
                    .arg("--db", &diamond_db),
            ])
            .name("databases"),
        );

    // ---------------------------------------------------------------- nail

    for &s in &args.nail_s {
        for (mode, query) in [("prf", &query_hmm), ("seq", &query_fa)] {
            let name = format!("nail-s{s:.1}-ms{MAX_SEQS}.{mode}");

            pl = pl.step(
                Step::serial([Cmd::new(&nail_bin)
                    .sub("search")
                    .arg("--mmseqs-path", &mmseqs_bin)
                    .arg("-t", args.threads)
                    .arg("--tmp-dir", tmp.join(&name))
                    .flag("--allow-overwrite")
                    .arg("--mmseqs-s", format!("{s:.1}"))
                    .arg("--mmseqs-max-seqs", MAX_SEQS)
                    .arg("-E", EVALUE)
                    .arg("--tbl-out", results.join(format!("{name}.tbl")))
                    .path(query)
                    .path(&target_fa)
                    .field("name", &name)])
                .name(name),
            );
        }
    }

    // -------------------------------------------------------------- hmmer

    for (mode, program, query) in [
        ("prf", &hmmsearch_bin, &query_hmm),
        ("seq", &phmmer_bin, &query_fa),
    ] {
        let name = format!("hmmer.{mode}");
        let kind = if mode == "prf" {
            bioio::split::Kind::Hmm
        } else {
            bioio::split::Kind::Fasta
        };

        let parts_dir = tmp.join(format!("{name}-parts"));
        let parts = bioio::split::write_splits(query, kind, hmmer_jobs, &parts_dir)?;

        pl = pl
            .step(
                Step::batched(
                    parts.len(),
                    parts.iter().enumerate().map(|(i, part)| {
                        Cmd::new(program)
                            .name(i.to_string())
                            .arg("--cpu", HMMER_CPU)
                            .arg("-E", EVALUE)
                            .arg("-o", "/dev/null")
                            .arg("--tblout", parts_dir.join(format!("{i}.tbl")))
                            .arg("--domtblout", parts_dir.join(format!("{i}.domtbl")))
                            .path(part)
                            .path(&target_fa)
                    }),
                )
                .name(&name),
            )
            .step(
                Step::serial([cat(
                    (0..parts.len()).map(|i| parts_dir.join(format!("{i}.tbl"))),
                    results.join(format!("{name}.tbl")),
                )
                .field("name", &name)])
                .name(format!("{name}.cat")),
            );
    }

    // ------------------------------------------------------------- mmseqs

    for &s in &args.mmseqs_s {
        for (mode, query_db) in [("prf", &query_prf_db), ("seq", &query_seq_db)] {
            let name = format!("mmseqs-s{s:.1}-ms{MAX_SEQS}.{mode}");
            let aln_db = tmp.join(format!("{name}/alnDB"));
            let work = tmp.join(format!("{name}/work"));

            pl = pl.step(
                Step::serial([
                    // mmseqs aborts rather than overwrite an existing
                    // alignment db, so a leftover from an earlier run has to
                    // go before this one can start
                    Cmd::new("rm")
                        .name("clean")
                        .flag("-rf")
                        .path(aln_db.parent().unwrap())
                        .path(&work),
                    Cmd::new("mkdir")
                        .name("dirs")
                        .flag("-p")
                        .path(aln_db.parent().unwrap())
                        .path(&work),
                    Cmd::new(&mmseqs_bin)
                        .name("search")
                        .sub("search")
                        .path(query_db)
                        .path(&target_db)
                        .path(&aln_db)
                        .path(&work)
                        .arg("--threads", args.threads)
                        .arg("-s", format!("{s:.1}"))
                        .arg("--max-seqs", MAX_SEQS)
                        .arg("-e", EVALUE),
                    Cmd::new(&mmseqs_bin)
                        .name("convertalis")
                        .sub("convertalis")
                        .path(query_db)
                        .path(&target_db)
                        .path(&aln_db)
                        .path(results.join(format!("{name}.tbl")))
                        .arg("--format-mode", 0)
                        .field("name", &name),
                ])
                .name(name),
            );
        }
    }

    // -------------------------------------------------------------- blast

    // sequence mode: one blastp call. deliberately no -evalue: matching the
    // other tools' 1e9 makes blast dramatically slower for no extra recall
    pl = pl.step(
        Step::serial([Cmd::new(&blastp_bin)
            .arg("-query", &query_fa)
            .arg("-db", &blast_db)
            .arg("-out", results.join("blast.seq.tbl"))
            .arg("-outfmt", 6)
            .arg("-num_threads", args.threads)
            .field("name", "blast.seq")])
        .name("blast.seq"),
    );

    // profile mode: psiblast takes one alignment at a time, so a run is one
    // invocation per family, output collected into a single table
    let blast_prf_tbl = results.join("blast.prf.tbl");
    pl = pl.step(
        Step::serial(afa.iter().enumerate().map(|(i, msa)| {
            let cmd = Cmd::new(&psiblast_bin)
                .name(
                    msa.file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("family")
                        .to_string(),
                )
                .arg("-in_msa", msa)
                .arg("-db", &blast_db)
                .arg("-outfmt", 6)
                .arg("-num_threads", args.threads)
                .arg("-comp_based_stats", 1)
                .arg("-num_iterations", 1)
                // one family's search is a fraction of the run, not a
                // separate run of its own, so its wall time is summed into
                // blast.prf's rather than reported as its own name
                .field("name", "blast.prf");

            // the first invocation truncates whatever an earlier run left
            // behind; the rest append to it
            if i == 0 {
                cmd.stdout_to(&blast_prf_tbl)
            } else {
                cmd.stdout(pail::Output::Append(blast_prf_tbl.clone()))
            }
        }))
        .name("blast.prf")
        .on_error(pail::OnError::Continue),
    );

    // ---------------------------------------------------------------- last

    pl = pl.step(
        Step::serial([Cmd::new(&lastal_bin)
            .path(&last_db)
            .path(&query_fa)
            .arg("-f", "BlastTab")
            .arg("-P", args.threads)
            .arg("-E", EVALUE)
            .stdout_to(results.join("last.seq.tbl"))
            .field("name", "last.seq")])
        .name("last.seq"),
    );

    // ------------------------------------------------------------- diamond

    for preset in ["default", "ultra-sensitive"] {
        let name = format!("diamond-{preset}.seq");
        let mut cmd = Cmd::new(&diamond_bin)
            .sub("blastp")
            .arg("--query", &query_fa)
            .arg("--db", &diamond_db)
            .arg("--out", results.join(format!("{name}.tbl")))
            .arg("--outfmt", 6)
            .arg("--threads", args.threads)
            .arg("--evalue", EVALUE)
            .field("name", &name);

        if preset == "ultra-sensitive" {
            cmd = cmd.flag("--ultra-sensitive");
        }

        pl = pl.step(Step::serial([cmd]).name(name));
    }

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(results.join("runs.tbl")))
        .build()
        .context("failed to build the run")?;

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
