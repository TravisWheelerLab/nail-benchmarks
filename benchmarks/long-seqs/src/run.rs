//! Searches each query against its paired target, one pair at a time.
//!
//! There is nothing to sweep and nothing to compare against: the question is
//! what fraction of the DP matrix nail computes at these lengths, and one run
//! of nail per pair answers it. So the pipeline is a step per pair, in order.
//!
//! Each search carries the fields `parse` reads back -- `name` for the stem of
//! its table, `tool` for how to read it, and `shard` for which pair it was --
//! the same contract the other two benchmarks write.

use std::path::PathBuf;

use anyhow::Context;
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step, Table};
use util::ledger;
use util::manifest;
use util::tools::nail;

use util::set::Set;

/// The name every pair's hit table and manifest row is filed under. One tool,
/// one setting, so every pair is the same run against a different target.
pub const RUN_NAME: &str = "nail";

/// Where the run's output lives, the same shape the other two benchmarks use.
pub struct Dirs {
    pub root: PathBuf,
    pub results: PathBuf,
    pub tmp: PathBuf,
}

impl Dirs {
    pub fn new(run: &std::path::Path, tmp: &std::path::Path) -> Dirs {
        Dirs {
            results: run.join("results"),
            tmp: tmp.to_owned(),
            root: run.to_owned(),
        }
    }

    pub fn table(&self, pair: &str) -> PathBuf {
        manifest::table_path(&self.results, RUN_NAME, pair)
    }
}

#[derive(Parser, Debug)]
pub struct Args {
    #[arg(short, long, default_value_t = 24)]
    pub threads: usize,

    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    #[arg(long)]
    pub tmp: Option<PathBuf>,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_as(&paths.set, &util::set::shape::PAIRS)?;

    let nail_bin = nail()?;

    let mut pl = PipelineBuilder::new().step(
        Cmd::new("mkdir")
            .name("dirs")
            .flag("-p")
            .path(&dirs.results),
    );

    for unit in set.units() {
        let pair = unit.name().to_string();

        pl = pl.step(
            Step::serial([Cmd::new(&nail_bin)
                .sub("search")
                .arg("-t", args.threads)
                .arg("--tmp-dir", dirs.tmp.join(&pair))
                .flag("--allow-overwrite")
                // widens the sparse band enough that these pairs align at all
                .arg("--f32-p", 5)
                .arg("--tbl-out", dirs.table(&pair))
                .path(unit.query_fa()?)
                .path(unit.target()?)
                .field(manifest::NAME, RUN_NAME)
                .field(manifest::TOOL, "nail")
                .field(manifest::SHARD, &pair)])
            .name(pair),
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
