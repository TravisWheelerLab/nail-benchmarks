//! One nail search per unit, at the settings every benchmark shares, with
//! `-s` on and stdout kept.
//!
//! michi discards a command's stdout unless told where to put it, and nail
//! prints its stage tree there, so the search names a file for it beside its
//! hit table.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step, Table};
use util::ledger;
use util::manifest;
use util::search::{Dirs, EVALUE, SEED_MODE, SEED_S, tag};
use util::set::Set;
use util::tools;

/// The name every unit's search is filed under: one tool, one setting.
pub const RUN_NAME: &str = "nail";

/// Where nail's stdout for one unit goes.
pub fn stats_path(dirs: &Dirs, unit: &str) -> PathBuf {
    dirs.table(RUN_NAME, unit).with_extension("stats")
}

#[derive(Parser, Debug)]
pub struct Args {
    /// Threads per search
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    #[arg(long)]
    pub tmp: Option<PathBuf>,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let threads = args.threads.context("--threads is required")?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_as(&paths.set, &util::set::shape::UNION)?;
    ensure!(set.units().count() > 0, "{} is empty", paths.set.display());

    let nail = tools::nail()?;
    let mmseqs = tools::mmseqs()?;

    let mut pl = PipelineBuilder::new().pool(threads).step(dirs.mkdir());

    for unit in set.units() {
        let name = unit.name().to_string();

        let cmd = Cmd::new(&nail)
            .sub("search")
            .arg("--mmseqs-path", &mmseqs)
            .arg("-t", threads)
            .arg("--tmp-dir", dirs.tmp.join(&name))
            .flag("--allow-overwrite")
            .arg("-E", EVALUE)
            .arg("--mmseqs-s", SEED_S)
            .arg("--seed-mode", SEED_MODE)
            .flag("-s")
            .arg("--tbl-out", dirs.table(RUN_NAME, &name))
            .path(unit.query_hmm()?)
            .path(unit.target()?)
            .stdout_to(stats_path(&dirs, &name));

        let cmd = tag(cmd, RUN_NAME, "nail", &[(manifest::SHARD, name.clone())]);
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

    // the ledger describes the results this run is about to replace, so it
    // goes before the run rather than after the failure of one
    ledger::clear(&dirs.root);
    pipeline.run()?;
    ledger::record(&dirs.root)
}
