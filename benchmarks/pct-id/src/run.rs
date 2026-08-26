//! Runs every tool against the one query/target set `pct-id build` made.
//!
//! nail and mmseqs sweep their prefilter sensitivity (`-s`); every other knob
//! is fixed. Every tool but last and diamond runs in both profile mode
//! (against `query.hmm`/`query.sto`) and sequence mode (against `query.fa`),
//! since the whole point of this benchmark is comparing the two.
//!
//! blast is deliberately left at its default E-value: raising it to match the
//! others makes it dramatically slower for no extra recall.

use std::path::PathBuf;

use anyhow::{Context, bail, ensure};
use clap::Parser;

use bench::manifest;
use bioio::split::Kind;
use pail::{Cmd, OnError, Output, PipelineBuilder, Progress, Step, Table};

use crate::inputs::Inputs;
use crate::search::{self, Bins, Dirs, MODE, PRF, SEQ, Split};

#[derive(Parser, Debug)]
pub struct Args {
    /// Which input set to search, naming inputs/<size>/.
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
    pub tmp: Option<PathBuf>,

    #[arg(long)]
    pub dry_run: bool,
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

    let set = Inputs::new(&args.size);
    if !set.exists() {
        bail!(
            "{} does not exist; run `pct-id build --size {}` first",
            set.dir().display(),
            args.size
        );
    }

    let mut dirs = Dirs::new(&set);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let (query_hmm, query_fa) = (set.query_hmm(), set.query_fa());
    let target_fa = set.target_fa();

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
                Cmd::new(&bins.mmseqs)
                    .name("createdb-target")
                    .sub("createdb")
                    .path(&target_fa)
                    .path(&target_db),
                Cmd::new(&bins.mmseqs)
                    .name("createdb-query-seq")
                    .sub("createdb")
                    .path(&query_fa)
                    .path(&query_seq_db),
                Cmd::new(&bins.mmseqs)
                    .name("convertmsa")
                    .sub("convertmsa")
                    .path(set.query_sto())
                    .path(&msa_db)
                    .arg("--identifier-field", 0),
                Cmd::new(&bins.mmseqs)
                    .name("msa2profile")
                    .sub("msa2profile")
                    .path(&msa_db)
                    .path(&query_prf_db)
                    .arg("--match-mode", 1),
                Cmd::new(&bins.makeblastdb)
                    .name("makeblastdb")
                    .arg("-in", &target_fa)
                    .arg("-dbtype", "prot")
                    .arg("-out", &blast_db),
                Cmd::new(&bins.lastdb)
                    .name("lastdb")
                    .arg("-p", &last_db)
                    .path(&target_fa),
                Cmd::new(&bins.diamond)
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
            let name = format!("nail-s{s:.1}-ms{}.{mode}", search::MAX_SEQS);

            pl = pl.step(
                Step::serial([Cmd::new(&bins.nail)
                    .sub("search")
                    .arg("--mmseqs-path", &bins.mmseqs)
                    .arg("-t", args.threads)
                    .arg("--tmp-dir", dirs.tmp.join(&name))
                    .flag("--allow-overwrite")
                    .arg("--mmseqs-s", format!("{s:.1}"))
                    .arg("--mmseqs-max-seqs", search::MAX_SEQS)
                    .arg("-E", search::EVALUE)
                    .arg("--tbl-out", dirs.table(&name))
                    .path(query)
                    .path(&target_fa)
                    .field(manifest::NAME, &name)
                    .field(manifest::TOOL, "nail")
                    .field(MODE, mode)
                    .field("s", format!("{s:.1}"))])
                .name(name),
            );
        }
    }

    // -------------------------------------------------------------- hmmer
    //
    // hmmsearch reads profiles and phmmer reads sequences, so the two modes are
    // two programs rather than one program with a flag.

    for (mode, tool, program, query, kind) in [
        (PRF, "hmmer", &bins.hmmsearch, &query_hmm, Kind::Hmm),
        (SEQ, "phmmer", &bins.phmmer, &query_fa, Kind::Fasta),
    ] {
        let name = format!("hmmer.{mode}");
        let split = Split::new(
            query,
            kind,
            dirs.tmp.join(format!("{name}-parts")),
            search::jobs(args.threads),
        );

        let hmmer = search::hmmer(
            program,
            &split,
            &dirs,
            search::Run {
                name: &name,
                tool,
                mode,
            },
            &target_fa,
            args.threads,
        );

        pl = pl
            .step(split.step(&name))
            .step(hmmer.search)
            .step(hmmer.cat);
    }

    // ------------------------------------------------------------- mmseqs

    for &s in &args.mmseqs_s {
        for (mode, query_db) in [(PRF, &query_prf_db), (SEQ, &query_seq_db)] {
            let name = format!("mmseqs-s{s:.1}-ms{}.{mode}", search::MAX_SEQS);

            let cmds = search::Mmseqs {
                bin: &bins.mmseqs,
                query_db,
                target_db: &target_db,
                scratch: dirs.tmp.join(&name),
                out: dirs.table(&name),
                threads: args.threads,
                s: format!("{s:.1}"),
            }
            .cmds();

            let searched = cmds.search.map(|cmd| {
                cmd.field(manifest::NAME, &name)
                    .field(manifest::TOOL, "mmseqs")
                    .field(MODE, mode)
                    .field("s", format!("{s:.1}"))
            });

            pl = pl.step(Step::serial(cmds.prep.into_iter().chain(searched)).name(name));
        }
    }

    // -------------------------------------------------------------- blast

    // sequence mode: one blastp call. deliberately no -evalue: matching the
    // other tools' 1e9 makes blast dramatically slower for no extra recall
    pl = pl.step(
        Step::serial([Cmd::new(&bins.blastp)
            .arg("-query", &query_fa)
            .arg("-db", &blast_db)
            .arg("-out", dirs.table("blast.seq"))
            .arg("-outfmt", 6)
            .arg("-num_threads", args.threads)
            .field(manifest::NAME, "blast.seq")
            .field(manifest::TOOL, "blast")
            .field(MODE, SEQ)])
        .name("blast.seq"),
    );

    // profile mode: psiblast takes one alignment at a time, so a run is one
    // invocation per family, output collected into a single table
    let blast_prf_tbl = dirs.table("blast.prf");
    pl = pl.step(
        Step::serial(set.afa_files()?.iter().enumerate().map(|(i, msa)| {
            let cmd = Cmd::new(&bins.psiblast)
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
                // one family's search is a fraction of the run, not a separate
                // run of its own, so its wall time is summed into blast.prf's
                // rather than reported as its own name
                .field(manifest::NAME, "blast.prf")
                .field(manifest::TOOL, "blast")
                .field(MODE, PRF);

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

    pl = pl.step(
        Step::serial([Cmd::new(&bins.lastal)
            .path(&last_db)
            .path(&query_fa)
            .arg("-f", "BlastTab")
            .arg("-P", args.threads)
            .arg("-E", search::EVALUE)
            .stdout_to(dirs.table("last.seq"))
            .field(manifest::NAME, "last.seq")
            .field(manifest::TOOL, "last")
            .field(MODE, SEQ)])
        .name("last.seq"),
    );

    // ------------------------------------------------------------- diamond

    for preset in ["default", "ultra-sensitive"] {
        let name = format!("diamond-{preset}.seq");
        let mut cmd = Cmd::new(&bins.diamond)
            .sub("blastp")
            .arg("--query", &query_fa)
            .arg("--db", &diamond_db)
            .arg("--out", dirs.table(&name))
            .arg("--outfmt", 6)
            .arg("--threads", args.threads)
            .arg("--evalue", search::EVALUE)
            .field(manifest::NAME, &name)
            .field(manifest::TOOL, "diamond")
            .field(MODE, SEQ)
            .field("preset", preset);

        if preset == "ultra-sensitive" {
            cmd = cmd.flag("--ultra-sensitive");
        }

        pl = pl.step(Step::serial([cmd]).name(name));
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
