//! Draws the pruning heatmaps, by handing summary.tbl to matplotlib.
//!
//! Only cloud-search has a surface to draw: the figure is a grid of -A against
//! -B, and a pipeline that swept something else has nothing to put in it. The
//! script says so and stops rather than failing.
//!
//! The drawing is python because matplotlib is what the other benchmarks plot
//! with and there is no reason for this one to be different. It goes through
//! the pipeline like everything else, which is what gets it a --dry-run, its
//! stderr kept on failure, and a line in the progress output.

use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, bail};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step};

const SCRIPT: &str = "scripts/plot.py";

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to draw. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Where the pdfs go. Defaults to what the label names
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,

    /// The interpreter to run the script with
    #[arg(long, default_value = "python3", value_name = "python")]
    python: String,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let summary = paths.analysis.join("summary.tbl");

    if !summary.is_file() {
        bail!(
            "no summary.tbl at {}; run `mgy parse summary` first",
            summary.display()
        );
    }

    let out = args.out.clone().unwrap_or_else(|| paths.figures.clone());

    let script = crate::dir().join(SCRIPT);
    if !script.is_file() {
        bail!("no plot script at {}", script.display());
    }

    // checked here rather than left to the pipeline, because a missing
    // matplotlib comes back as a python traceback in a stderr file rather than
    // as anything that reads like an answer
    matplotlib(&args.python)?;

    // the script goes in the subcommand slot rather than in a path: a Cmd
    // renders its options before its positionals, and python needs the script
    // ahead of everything
    let cmd = Cmd::new(&args.python)
        .name("plot")
        .sub(script.to_string_lossy())
        .arg("--out", &out)
        .path(&summary);

    let pipeline = PipelineBuilder::new()
        .step(Step::serial([cmd]))
        .stderr_dir(paths.tmp.join("plot-stderr"))
        .sink(Progress::new())
        .build()?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()?;

    println!("\nfigures in {}", out.display());
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
