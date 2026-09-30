//! Which stage of the nail pipeline drops the hits hmmer finds.
//!
//! Runs hmmer once, then one nail per sensitivity, each seeded static with no
//! cap so it aligns everything its prefilter returned. There is no grid here
//! the way there is in cloud-search: an arm is the most nail could get at that
//! sensitivity, and the question is where the hits hmmer found fall out of it.
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

use util::ledger;
use util::manifest;
use util::search::{Bins, Dirs, SEED_S, Split};
use util::set::Set;
use util::split::Kind;

/// The column hmmer's run becomes, which every arm is measured against.
const HMMER: &str = "hmmer";

/// One arm of the sweep: a sensitivity, and the run name its column takes.
//
// static seeding with no cap, so the arm aligns everything the prefilter
// returned at that sensitivity and the seed list is the most nail could get
// from it. hmmer and the query split sit outside: the truth set is the same
// for every arm, and it is most of the wall clock
struct Arm {
    name: String,
    s: String,
    /// `--prog-n` and `--prog-f` for a prog arm, which seeds at the shared
    /// sensitivity and stops per query; none for a static arm.
    prog: Option<(usize, f64)>,
}

impl Arm {
    fn at(s: f64) -> Arm {
        Arm {
            name: format!("s{s}"),
            s: format!("{s}"),
            prog: None,
        }
    }

    /// The point the static ceiling is read against: what nail does by
    /// default, at the sensitivity every benchmark here seeds at.
    fn prog(n: usize, f: f64) -> Arm {
        Arm {
            name: format!("prog-n{n}-f{f}"),
            s: SEED_S.to_string(),
            prog: Some((n, f)),
        }
    }
}

/// nail's own unbounded `--mmseqs-max-seqs`, the default it uses in prog mode.
const UNBOUNDED: u32 = i32::MAX as u32;

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// nail's -E, set far above its default so the final e-value gate can't be
    /// mistaken for a cloud/align filter. Only lower this to study the e-value
    /// gate itself
    #[arg(long, default_value_t = 1e6, value_name = "X")]
    nail_evalue: f64,

    /// `--mmseqs-s` values to sweep, each seeded static and unbounded
    #[arg(
        long = "s",
        value_delimiter = ',',
        default_value = "12.0,10.0,7.5",
        value_name = "X,X,..."
    )]
    sensitivities: Vec<f64>,

    /// A prog arm to add, as `n:f`, seeding at the shared sensitivity
    #[arg(long, value_delimiter = ',', value_name = "N:F,...")]
    prog: Vec<String>,

    /// Threads per search, and the cores each search is pinned to
    #[arg(short, long)]
    threads: Option<usize>,

    #[arg(long)]
    tmp: Option<PathBuf>,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let threads = args.threads.context("--threads is required")?;

    ensure!(
        threads.is_multiple_of(util::search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        util::search::HMMER_CPU
    );

    let bins = Bins::find()?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_needing(&paths.set, crate::NEEDS)?;
    ensure!(set.units().count() > 0, "{} is empty", paths.set.display());

    // nail creates the last component of its --tmp-dir and no more, so the
    // directory the arms sit under has to exist before the first one runs
    let mut pl = PipelineBuilder::new()
        .pool(threads)
        .step(dirs.mkdir().path(dirs.tmp.join("align")));

    let mut arms: Vec<Arm> = args.sensitivities.iter().map(|&s| Arm::at(s)).collect();
    for spec in &args.prog {
        let (n, f) = spec
            .split_once(':')
            .with_context(|| format!("--prog wants n:f, got {spec:?}"))?;
        arms.push(Arm::prog(n.parse()?, f.parse()?));
    }

    ensure!(!arms.is_empty(), "the sweep has no arms");
    println!(
        "{} arms over {} {}",
        arms.len(),
        set.units().count(),
        if set.units().count() == 1 {
            "unit"
        } else {
            "units"
        }
    );

    // every unit gets every arm. On a `fixed` set that is the one shard the
    // label names; on a `cross` it is each target source in turn
    for unit in set.units() {
        let shard = unit.name().to_string();
        let query_hmm = unit.query_hmm()?;
        let target = unit.target()?;

        // per unit rather than per run: a cross can pair more than one query
        // source, and two of them cut into the same directory would search each
        // other's parts
        let split = Split::new(
            &query_hmm,
            Kind::Hmm,
            dirs.tmp.join("hmmer-query").join(&shard),
            util::search::jobs(threads),
        );
        pl = pl.step(split.step("split", &[(manifest::SHARD, shard.clone())]));

        let hmmer = util::search::hmmer(
            &bins.hmmsearch,
            &split,
            &dirs,
            HMMER,
            "hmmer",
            &shard,
            &target,
            util::search::EVALUE,
            true,
            &[],
        );
        pl = pl.step(hmmer.search).step(hmmer.cat);

        // one nail per arm: static seeding with no cap aligns everything
        // the prefilter returned, so there is no list to replay. the seed
        // list is written beside the table for parse to read which pairs
        // the seeding offered, and mmseqs' databases stay in the tmp dir
        // for depth to read a pair's rank
        for arm in &arms {
            let mut cmd = Cmd::new(&bins.nail)
                .sub("search")
                .arg("--mmseqs-path", &bins.mmseqs)
                .arg("-t", threads)
                .arg(
                    "--tmp-dir",
                    crate::scores::depth::prefilter_dir(&dirs.root, &arm.name, &shard),
                )
                .arg("--mmseqs-s", &arm.s)
                .arg("--mmseqs-max-seqs", UNBOUNDED);

            cmd = match arm.prog {
                Some((n, f)) => cmd
                    .arg("--seed-mode", "prog")
                    .arg("--prog-n", n)
                    .arg("--prog-f", f),
                None => cmd.arg("--seed-mode", "static"),
            };

            let cmd = cmd
                .arg("--seeds-out", dirs.seeds(&arm.name, &shard))
                .arg("-E", args.nail_evalue)
                .arg("--tbl-out", dirs.table(&arm.name, &shard))
                .flag("--allow-overwrite")
                .path(&query_hmm)
                .path(&target)
                .field(manifest::NAME, &arm.name)
                .field(manifest::TOOL, "nail")
                // the list this run wrote, which is what it aligned
                .field(manifest::SEEDS, &arm.name)
                .field(manifest::SHARD, &shard)
                .field("s", &arm.s)
                .field("E", args.nail_evalue);

            pl = pl.step(
                Step::serial([cmd])
                    .name(format!("{}.{shard}", arm.name))
                    .cores(threads),
            );
        }
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
