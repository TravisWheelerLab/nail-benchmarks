//! Draws `cells.pdf` by handing every form's table to matplotlib, through the
//! pipeline so a failed script keeps its stderr.

use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, bail};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step};
use util::set::{Set, shape};

use crate::parse::table_name;

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Where the pdf goes. Defaults to what the label names
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,

    /// The interpreter to run the script with
    #[arg(long, default_value = "python3", value_name = "python")]
    python: String,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let set = Set::load_as(&paths.set, &shape::FORMS)?;
    let dir = args.out.clone().unwrap_or_else(|| paths.figures.clone());

    let script = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scripts/plot.py");

    // the script reads its tables in the order the forms are given, so a
    // form's colour follows its name rather than its row in the set
    let mut cmd = Cmd::new(&args.python)
        .name("cells")
        .sub(script.to_string_lossy());
    for unit in set.units() {
        let form = unit.need("form")?;
        let table = paths.analysis.join(table_name(form));
        if !table.is_file() {
            bail!("no {}; run `cells parse` first", table.display());
        }
        cmd = cmd.arg(format!("--{form}"), &table);
    }
    let cmd = cmd.path(dir.join("cells.pdf"));

    // checked here rather than left to the pipeline, because a missing
    // matplotlib comes back as a python traceback in a stderr file
    matplotlib(&args.python)?;
    std::fs::create_dir_all(&dir)?;

    let pipeline = PipelineBuilder::new()
        .step(Step::serial([cmd]).name("cells"))
        .stderr_dir(paths.tmp.join("plot-stderr"))
        .sink(Progress::new())
        .build()
        .context("failed to build the figure")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()?;

    println!("\nfigure in {}", dir.display());
    Ok(())
}

/// Checks that the interpreter is there and can import matplotlib.
fn matplotlib(python: &str) -> anyhow::Result<()> {
    let out = Command::new(python)
        .args(["-c", "import matplotlib"])
        .output()
        .with_context(|| format!("couldn't run {python}"))?;

    if !out.status.success() {
        bail!(
            "{python} can't import matplotlib:\n{}",
            String::from_utf8_lossy(&out.stderr).trim()
        );
    }

    Ok(())
}
