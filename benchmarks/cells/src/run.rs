//! One nail search per form, at nail's defaults with the reporting
//! threshold opened.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{PipelineBuilder, Progress, Table};
use util::ledger;
use util::manifest;
use util::search::{self, Dirs};
use util::set::{Set, shape};
use util::tools;

/// The name every form's search is filed under.
pub const RUN_NAME: &str = "nail";

/// The E-value nail reports down to: high enough that every pair through
/// nail's filters is reported, since a pair with no signal is the point of
/// two of the forms.
const EVALUE: &str = "1e9";

/// nail's arguments beyond `-E`. `--f32-p 5` prints the cell fraction to five
/// decimals rather than three, which would round the smallest to 0.
const ARGS: [(&str, &str); 1] = [("--f32-p", "5")];

#[derive(Parser, Debug)]
pub struct Args {
    /// Threads per search, and the cores the run is pinned to
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Where nail's scratch goes. Defaults to what the label names
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

    let set = Set::load_as(&paths.set, &shape::FORMS)?;
    ensure!(set.units().count() > 0, "{} is empty", paths.set.display());

    let nail = tools::nail()?;
    let mmseqs = tools::mmseqs()?;

    let mut pl = PipelineBuilder::new().pool(threads).step(dirs.mkdir());

    for unit in set.units() {
        let name = unit.name();

        let step = search::nail(
            &nail,
            &mmseqs,
            &unit.query_hmm()?,
            &unit.target()?,
            &dirs.table(RUN_NAME, name),
            &dirs.tmp.join(name),
            threads,
            EVALUE,
            &ARGS.map(|(flag, value)| (flag, value.to_string())),
            &[],
            RUN_NAME,
            &[(manifest::SHARD, name.to_string())],
        );
        pl = pl.step(step.name(format!("{RUN_NAME}.{name}")));
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
