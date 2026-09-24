//! Runs every tool against the one query/target set `build-set` made.
//!
//! nail and mmseqs sweep their prefilter sensitivity (`-s`); every other knob
//! is fixed. Every tool but last and diamond runs in both profile mode
//! (against `query.hmm`/`query.sto`) and sequence mode (against `query.fa`),
//! since the whole point of this benchmark is comparing the two.
//!
//! blast is deliberately left at its default E-value: raising it to match the
//! others makes it dramatically slower for no extra recall.

use std::path::{Path, PathBuf};

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Cmd, OnError, Output, PipelineBuilder, Progress, Step, Table};
use util::ledger;
use util::split::Kind;

use util::search::{Dirs, Mmseqs, Split, hmmer, jobs, nail, tag};

/// The `mode` field, and its two values: a profile searched a target, or a
/// sequence did.
pub const MODE: &str = "mode";
pub const PRF: &str = "prf";
pub const SEQ: &str = "seq";

/// Every tool reports down to here. A benchmark drawing an ROC curve wants
/// every hit, not the ones that clear a threshold.
pub const EVALUE: &str = "1e9";

/// How many targets the prefilters promote. mmseqs' own default is 300, which
/// loses hits nail's seeding keeps.
const MAX_SEQS: usize = 2000;

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// nail's --mmseqs-s values to sweep
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "5.7,7.5,10.0,12.0,14.0",
        value_name = "X,X,..."
    )]
    pub nail_s: Vec<f32>,

    /// mmseqs' -s values to sweep
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
    pub tmp: Option<PathBuf>,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(util::search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        util::search::HMMER_CPU
    );
    ensure!(!args.nail_s.is_empty(), "--nail-s needs at least one value");
    ensure!(
        !args.mmseqs_s.is_empty(),
        "--mmseqs-s needs at least one value"
    );

    // resolved before anything runs: a sweep that dies two tools in has spent
    // its wall time on numbers that are no longer comparable to the ones it
    // did not reach
    let nail_bin = util::tools::nail()?;
    let mmseqs = util::tools::mmseqs()?;
    let hmmsearch = util::tools::hmmsearch()?;
    let phmmer = util::tools::phmmer()?;
    let blastp_bin = util::tools::blastp()?;
    let psiblast_bin = util::tools::psiblast()?;
    let makeblastdb = util::tools::makeblastdb()?;
    let lastal_bin = util::tools::lastal()?;
    let lastdb = util::tools::lastdb()?;
    let diamond_bin = util::tools::diamond()?;

    // Set::load_as is the check: it refuses a set that is not there and one
    // that is the wrong shape, before any tool runs
    let inp = crate::Inputs::open(&paths.set)?;

    let mut dirs = Dirs::new(paths.search(), paths.tmp.join("search"));
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let (query_hmm, query_fa) = (inp.query_hmm.clone(), inp.query_fa.clone());
    let target_fa = inp.target_fa.clone();

    let mmseqs_dir = dirs.tmp.join("mmseqs");
    let target_db = mmseqs_dir.join("targetDB/targetDB");
    let query_prf_db = mmseqs_dir.join("queryDB-prf/queryDB");
    let query_seq_db = mmseqs_dir.join("queryDB-seq/queryDB");
    let msa_db = mmseqs_dir.join("msaDB/msaDB");
    let blast_db = dirs.tmp.join("blast/db");
    let last_db = dirs.tmp.join("last/db");
    let diamond_db = dirs.tmp.join("diamond/db");

    let mut pl = PipelineBuilder::new()
        .step(dirs.clean())
        // a rejection settles a pair of the search it was drawn from, so a new
        // search leaves none behind to be read against it
        .step(
            Cmd::new("rm")
                .name("clean-reject")
                .flag("-rf")
                .path(paths.reject()),
        )
        .step(
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(target_db.parent().expect("targetDB has a parent"))
                .path(query_prf_db.parent().expect("queryDB has a parent"))
                .path(query_seq_db.parent().expect("queryDB has a parent"))
                .path(msa_db.parent().expect("msaDB has a parent"))
                .path(blast_db.parent().expect("blast db has a parent"))
                .path(last_db.parent().expect("last db has a parent"))
                .path(diamond_db.parent().expect("diamond db has a parent")),
        )
        .step(
            Step::serial([
                Cmd::new(&mmseqs)
                    .name("createdb-target")
                    .sub("createdb")
                    .path(&target_fa)
                    .path(&target_db),
                Cmd::new(&mmseqs)
                    .name("createdb-query-seq")
                    .sub("createdb")
                    .path(&query_fa)
                    .path(&query_seq_db),
                Cmd::new(&mmseqs)
                    .name("convertmsa")
                    .sub("convertmsa")
                    .path(&inp.query_sto)
                    .path(&msa_db)
                    .arg("--identifier-field", 0),
                Cmd::new(&mmseqs)
                    .name("msa2profile")
                    .sub("msa2profile")
                    .path(&msa_db)
                    .path(&query_prf_db)
                    .arg("--match-mode", 1),
                Cmd::new(&makeblastdb)
                    .name("makeblastdb")
                    .arg("-in", &target_fa)
                    .arg("-dbtype", "prot")
                    .arg("-out", &blast_db),
                Cmd::new(&lastdb)
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
        for (mode, query) in [(PRF, &query_hmm), (SEQ, &query_fa)] {
            let name = format!("nail-s{s:.1}-ms{}.{mode}", MAX_SEQS);

            pl = pl.step(nail(
                &nail_bin,
                &mmseqs,
                query,
                &target_fa,
                &dirs.table(&name, ""),
                &dirs.tmp.join(&name),
                args.threads,
                EVALUE,
                &[
                    ("--mmseqs-s", format!("{s:.1}")),
                    ("--mmseqs-max-seqs", MAX_SEQS.to_string()),
                ],
                &[],
                &name,
                &[(MODE, mode.to_string()), ("s", format!("{s:.1}"))],
            ));
        }
    }

    // -------------------------------------------------------------- hmmer
    //
    // hmmsearch reads profiles and phmmer reads sequences, so the two modes are
    // two programs rather than one program with a flag.

    for (mode, tool, program, query, kind) in [
        (PRF, "hmmer", &hmmsearch, &query_hmm, Kind::Hmm),
        (SEQ, "phmmer", &phmmer, &query_fa, Kind::Fasta),
    ] {
        let name = format!("hmmer.{mode}");
        let split = Split::new(
            query,
            kind,
            dirs.tmp.join(format!("{name}-parts")),
            jobs(args.threads),
        );

        let hmmer = hmmer(
            program,
            &split,
            &dirs,
            &name,
            tool,
            "",
            &target_fa,
            EVALUE,
            false,
            &[(MODE, mode.to_string())],
        );

        pl = pl
            .step(split.step(&name, &[]))
            .step(hmmer.search)
            .step(hmmer.cat);
    }

    // ------------------------------------------------------------- mmseqs

    for &s in &args.mmseqs_s {
        for (mode, query_db) in [(PRF, &query_prf_db), (SEQ, &query_seq_db)] {
            let name = format!("mmseqs-s{s:.1}-ms{}.{mode}", MAX_SEQS);

            let scratch = dirs.tmp.join(&name);

            let cmds = Mmseqs {
                bin: &mmseqs,
                query_db,
                target_db: &target_db,
                aln_db: scratch.join("alnDB"),
                work: scratch.join("work"),
                out: dirs.table(&name, ""),
                threads: args.threads,
                s: Some(format!("{s:.1}")),
                max_seqs: Some(MAX_SEQS),
                evalue: EVALUE,
            }
            .cmds();

            let prep = [
                Cmd::new("rm").name("clean").flag("-rf").path(&scratch),
                Cmd::new("mkdir").name("dirs").flag("-p").path(&scratch),
            ];
            let search = [cmds.search, cmds.convert];

            let fields = [(MODE, mode.to_string()), ("s", format!("{s:.1}"))];
            let searched = search.map(|cmd| tag(cmd, &name, "mmseqs", &fields));

            pl = pl.step(Step::serial(prep.into_iter().chain(searched)).name(name));
        }
    }

    // -------------------------------------------------------------- blast

    // sequence mode: one blastp call. deliberately no -evalue: matching the
    // other tools' 1e9 makes blast dramatically slower for no extra recall
    pl = pl.step(blastp(
        &blastp_bin,
        &query_fa,
        &blast_db,
        &dirs.table("blast.seq", ""),
        args.threads,
        "blast.seq",
    ));

    // profile mode: psiblast takes one alignment at a time, so a run is one
    // invocation per family, output collected into a single table
    let blast_prf_tbl = dirs.table("blast.prf", "");

    let afa = alignments(&inp.afa)?;

    pl = pl.step(
        Step::serial(afa.iter().enumerate().map(|(i, msa)| {
            // one family's search is a fraction of the run rather than a run
            // of its own, so every command carries blast.prf's name and their
            // wall times sum into it
            let cmd = psiblast(&psiblast_bin, msa, &blast_db, args.threads, "blast.prf");

            // the first invocation truncates whatever an earlier run left
            // behind; the rest append to it
            match i {
                0 => cmd.stdout_to(&blast_prf_tbl),
                _ => cmd.stdout(Output::Append(blast_prf_tbl.clone())),
            }
        }))
        .name("blast.prf")
        .on_error(OnError::Continue),
    );

    // ---------------------------------------------------------------- last

    pl = pl.step(lastal(
        &lastal_bin,
        &last_db,
        &query_fa,
        &dirs.table("last.seq", ""),
        args.threads,
        "last.seq",
    ));

    // ------------------------------------------------------------- diamond

    for preset in ["default", "ultra-sensitive"] {
        let name = format!("diamond-{preset}.seq");
        pl = pl.step(diamond(
            &diamond_bin,
            &query_fa,
            &diamond_db,
            &dirs.table(&name, ""),
            args.threads,
            &name,
            preset,
        ));
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

    // the ledger describes the results this run is about to replace, so it
    // goes before the run rather than after the failure of one
    ledger::clear(&dirs.root);
    pipeline.run()?;
    ledger::record(&dirs.root)
}

/// One aligned fasta per family, the input psiblast reads.
fn alignments(dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let mut afa: Vec<PathBuf> = std::fs::read_dir(dir)
        .with_context(|| format!("failed to read {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "afa"))
        .collect();

    // sorted so two runs search the families in the same order, which is what
    // makes their wall times comparable
    afa.sort();
    Ok(afa)
}

/// blastp over sequences: one call, its own table.
fn blastp(bin: &Path, query: &Path, db: &Path, out: &Path, threads: usize, name: &str) -> Step {
    // deliberately no -evalue: matching the other tools' 1e9 makes blast
    // dramatically slower for no extra recall
    let cmd = Cmd::new(bin)
        .arg("-query", query)
        .arg("-db", db)
        .arg("-out", out)
        .arg("-outfmt", 6)
        .arg("-num_threads", threads);

    Step::serial([tag(cmd, name, "blast", &[(MODE, SEQ.to_string())])]).name(name)
}

/// psiblast over one family's alignment. It takes an alignment at a time, so
/// a run is one of these per family and the caller collects them.
fn psiblast(bin: &Path, msa: &Path, db: &Path, threads: usize, name: &str) -> Cmd {
    let cmd = Cmd::new(bin)
        .name(
            msa.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("family")
                .to_string(),
        )
        .arg("-in_msa", msa)
        .arg("-db", db)
        .arg("-outfmt", 6)
        .arg("-num_threads", threads)
        .arg("-comp_based_stats", 1)
        .arg("-num_iterations", 1);

    // one family's search is a fraction of the run rather than a run of its
    // own, so every command carries the run's name and their wall times sum
    tag(cmd, name, "blast", &[(MODE, PRF.to_string())])
}

/// lastal over sequences, writing blast's tabular format to stdout.
fn lastal(bin: &Path, db: &Path, query: &Path, out: &Path, threads: usize, name: &str) -> Step {
    let cmd = Cmd::new(bin)
        .path(db)
        .path(query)
        .arg("-f", "BlastTab")
        .arg("-P", threads)
        .arg("-E", EVALUE)
        .stdout_to(out);

    Step::serial([tag(cmd, name, "last", &[(MODE, SEQ.to_string())])]).name(name)
}

/// diamond blastp at one of its presets.
fn diamond(
    bin: &Path,
    query: &Path,
    db: &Path,
    out: &Path,
    threads: usize,
    name: &str,
    preset: &str,
) -> Step {
    let mut cmd = Cmd::new(bin)
        .sub("blastp")
        .arg("--query", query)
        .arg("--db", db)
        .arg("--out", out)
        .arg("--outfmt", 6)
        .arg("--threads", threads)
        .arg("--evalue", EVALUE);

    if preset == "ultra-sensitive" {
        cmd = cmd.flag("--ultra-sensitive");
    }

    let fields = [(MODE, SEQ.to_string()), ("preset", preset.to_string())];
    Step::serial([tag(cmd, name, "diamond", &fields)]).name(name)
}
