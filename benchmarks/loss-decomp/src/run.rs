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

use michi::{Cmd, PipelineBuilder, Progress, Step, Table};

use search::sweeps::SEED_S;
use search::{self, Bins, Dirs, Split};
use util::ledger;
use util::manifest;
use util::set::Set;

/// The column hmmer's run becomes, which every arm is measured against.
const HMMER: &str = "hmmer";

/// One arm of the sweep: a seeding, and the run name its column takes.
//
// the knobs are all seeding knobs, so an arm is a seed list of its own and a
// nail that replays it. hmmer and the query split sit outside: the truth set
// is the same for every arm, and it is most of the wall clock
struct Arm {
    name: String,
    seeding: search::Seeding<'static>,
}

impl Arm {
    fn static_(max_seqs: usize) -> Arm {
        Arm {
            name: format!("static-ms{max_seqs}"),
            seeding: search::Seeding {
                mmseqs_s: SEED_S,
                mode: "static",
                max_seqs: Some(max_seqs),
                prog_n: None,
                prog_f: None,
            },
        }
    }

    fn prog(n: usize, f: f64) -> Arm {
        Arm {
            name: format!("prog-n{n}-f{f}"),
            seeding: search::Seeding {
                mmseqs_s: SEED_S,
                mode: "prog",
                max_seqs: None,
                prog_n: Some(n),
                prog_f: Some(f),
            },
        }
    }
}

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Which unit of the set to search
    #[arg(long, default_value = "1", value_name = "N")]
    shard: String,

    /// nail's -E, set far above its default so the final e-value gate can't be
    /// mistaken for a cloud/align filter. Only lower this to study the e-value
    /// gate itself
    #[arg(long, default_value_t = 1e6, value_name = "X")]
    nail_evalue: f64,

    /// `--mmseqs-max-seqs` values to sweep in static mode
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "200,2000,20000",
        value_name = "N,N,..."
    )]
    max_seqs: Vec<usize>,

    /// `--prog-n` values to sweep in prog mode
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "50,200,800",
        value_name = "N,N,..."
    )]
    prog_n: Vec<usize>,

    /// `--prog-f` values to sweep in prog mode
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "0.001,0.01,0.1",
        value_name = "X,X,..."
    )]
    prog_f: Vec<f64>,

    /// Threads per search, and the cores each search is pinned to
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    #[arg(long)]
    tmp: Option<PathBuf>,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        search::HMMER_CPU
    );

    let bins = Bins::find()?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_as(&paths.set, crate::SHAPE)?;

    let unit = set
        .units()
        .find(|u| u.name() == args.shard)
        .with_context(|| {
            format!(
                "{} has no unit {:?}; it holds {}",
                paths.set.display(),
                args.shard,
                set.units().map(|u| u.name()).collect::<Vec<_>>().join(", ")
            )
        })?;

    let query_hmm = unit.query_hmm()?;
    let target = unit.target()?;

    let split = Split::new(
        &query_hmm,
        dirs.tmp.join("hmmer-query"),
        search::jobs(args.threads),
    );

    let mut pl = PipelineBuilder::new()
        .step(dirs.mkdir())
        .step(split.step(&[]));

    // static aligns the whole prefilter, so max-seqs bounds it; prog aligns
    // from prog-n upward while the hit fraction holds. one arm, one seed list
    let arms: Vec<Arm> = args
        .max_seqs
        .iter()
        .map(|&n| Arm::static_(n))
        .chain(
            args.prog_n
                .iter()
                .flat_map(|&n| args.prog_f.iter().map(move |&f| Arm::prog(n, f))),
        )
        .collect();

    ensure!(!arms.is_empty(), "the sweep has no arms");
    println!("{} arms over shard {}", arms.len(), args.shard);

    for arm in &arms {
        pl = pl.step(
            search::seed(
                &bins.nail,
                &bins.mmseqs,
                &query_hmm,
                &target,
                &args.shard,
                &dirs.seeds(&arm.name),
                &dirs,
                args.threads,
                &arm.seeding,
                &[(manifest::STAGE, search::SEED.to_string())],
            )
            .name(format!("seeds.{}", arm.name)),
        );
    }

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

    for arm in &arms {
        pl = pl.step(
            Step::serial([Cmd::new(&bins.nail)
                .sub("search")
                // nail looks for mmseqs at startup even when it is replaying
                // seeds and will never call it, and nothing here is on PATH
                .arg("--mmseqs-path", &bins.mmseqs)
                .arg("-t", args.threads)
                .arg("--seeds", dirs.seeds(&arm.name))
                .arg("-E", args.nail_evalue)
                .arg("--tmp-dir", dirs.tmp.join("align").join(&arm.name))
                .arg("--tbl-out", dirs.table(&arm.name, &args.shard))
                .flag("--allow-overwrite")
                .path(&query_hmm)
                .path(&target)
                .field(manifest::NAME, &arm.name)
                .field(manifest::TOOL, "nail")
                .field(manifest::SHARD, &args.shard)
                .field("E", args.nail_evalue)])
            .name(arm.name.clone())
            .cores(args.threads),
        );
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
