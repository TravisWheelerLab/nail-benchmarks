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
use util::search::{Bins, Dirs, Split};
use util::set::Set;
use util::split::Kind;

/// The column hmmer's run becomes, which every arm is measured against.
const HMMER: &str = "hmmer";

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

    /// Sensitivities to add a full-DP arm at: `s<X>-full` replays `s<X>`'s
    /// seed list with `--full-dp`, so the only thing moving is the cloud
    #[arg(long = "full-dp", value_delimiter = ',', value_name = "X,X,...")]
    full_dp: Vec<f64>,

    /// Run only these arms, and add them to the ledger of a run that has
    /// finished rather than starting the record over. `hmmer` names hmmer
    #[arg(long, value_delimiter = ',', value_name = "name,...")]
    only: Vec<String>,

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

    // an arm is its name, its sensitivity, and whether it replays the
    // static arm's seed list under --full-dp
    let mut arms: Vec<(String, String, bool)> = args
        .sensitivities
        .iter()
        .map(|s| (format!("s{s}"), format!("{s}"), false))
        .collect();
    arms.extend(
        args.full_dp
            .iter()
            .map(|s| (format!("s{s}-full"), format!("{s}"), true)),
    );

    ensure!(!arms.is_empty(), "the sweep has no arms");

    let wanted = |name: &str| args.only.is_empty() || args.only.iter().any(|o| o == name);
    for name in &args.only {
        ensure!(
            name == HMMER || arms.iter().any(|(arm, _, _)| arm == name),
            "--only names {name:?}, which is not an arm of this run"
        );
    }
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
        if wanted(HMMER) {
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
        }

        // one nail per arm: static seeding with no cap aligns everything
        // the prefilter returned, so there is no list to replay. the seed
        // list is written beside the table for parse to read which pairs
        // the seeding offered, and mmseqs' databases stay in the tmp dir
        // for depth to read a pair's rank
        for (name, s, full) in &arms {
            if !wanted(name) {
                continue;
            }

            let cmd = Cmd::new(&bins.nail)
                .sub("search")
                // nail looks for mmseqs at startup even when it replays
                // seeds and never calls it
                .arg("--mmseqs-path", &bins.mmseqs)
                .arg("-t", threads);

            // the seeding it replays is the static arm at the same -s; a
            // static arm seeds itself and writes the list the full arm reads
            let seeding = format!("s{s}");
            let cmd = match full {
                true => cmd
                    .arg("--seeds", dirs.seeds(&seeding, &shard))
                    .flag("--full-dp")
                    .arg("--tmp-dir", dirs.tmp.join("full"))
                    .field("dp", "full"),
                false => cmd
                    .arg(
                        "--tmp-dir",
                        crate::scores::depth::prefilter_dir(&dirs.root, name, &shard),
                    )
                    .arg("--mmseqs-s", s)
                    .arg("--mmseqs-max-seqs", UNBOUNDED)
                    .arg("--seed-mode", "static")
                    .arg("--seeds-out", dirs.seeds(name, &shard)),
            };

            let cmd = cmd
                .arg("-E", args.nail_evalue)
                .arg("--tbl-out", dirs.table(name, &shard))
                .flag("--allow-overwrite")
                .path(&query_hmm)
                .path(&target)
                .field(manifest::NAME, name)
                .field(manifest::TOOL, "nail")
                // the list it aligned: its own, or the static arm's
                .field(manifest::SEEDS, &seeding)
                .field(manifest::SHARD, &shard)
                .field("s", s)
                .field("E", args.nail_evalue);

            pl = pl.step(
                Step::serial([cmd])
                    .name(format!("{name}.{shard}"))
                    .cores(threads),
            );
        }
    }

    // an added arm keeps its own manifest beside the run's, and joins the
    // run's ledger rather than replacing it
    let manifest = match args.only.is_empty() {
        true => dirs.root.join("manifest.tbl"),
        false => dirs
            .root
            .join(format!("manifest-{}.tbl", args.only.join("+"))),
    };

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(&manifest))
        .build()
        .context("failed to build the run")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    if !args.only.is_empty() {
        pipeline.run()?;
        return ledger::add(&dirs.root, &manifest);
    }

    // the ledger describes the results this run is about to replace, so it
    // goes before the run rather than after the failure of one
    ledger::clear(&dirs.root);
    pipeline.run()?;
    ledger::record(&dirs.root)
}
