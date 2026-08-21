//! Draws the figures, by handing grid.tbl to matplotlib.
//!
//! The drawing is python because matplotlib is what the other benchmarks plot
//! with and there is no reason for this one to be different. It goes through
//! the pipeline like everything else, which is what gets it a --dry-run, its
//! stderr kept on failure, and a line in the progress output.

use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, bail};
use clap::Parser;

use pipeline::{Cmd, PipelineBuilder, Progress, Step};

const SCRIPT: &str = "scripts/plot.py";

#[derive(Parser, Debug)]
pub struct Args {
    /// The grid.tbl `parse grid` wrote, or the benchmark directory holding one.
    grid: PathBuf,

    /// Where the pdfs go. Defaults to figures/ beside grid.tbl.
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,

    /// Draw only this figure. Repeatable; every one by default.
    #[arg(long, value_name = "NAME")]
    only: Vec<String>,

    /// The interpreter to run the script with.
    #[arg(long, default_value = "python3", value_name = "python")]
    python: String,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let grid = match args.grid.is_dir() {
        true => args.grid.join("grid.tbl"),
        false => args.grid.clone(),
    };

    if !grid.is_file() {
        bail!(
            "no grid.tbl at {}; run `cloud-search parse grid` first",
            grid.display()
        );
    }

    let out = match args.out {
        Some(dir) => dir,
        None => grid
            .parent()
            .context("grid.tbl has no directory")?
            .join("figures"),
    };

    let script = crate::dir().join(SCRIPT);
    if !script.is_file() {
        bail!("no plot script at {}", script.display());
    }

    // checked here rather than left to the pipeline, because a missing
    // matplotlib comes back as a python traceback in a stderr file rather than
    // as anything that reads like an answer
    matplotlib(&args.python)?;

    // the script goes in the subcommand slot rather than in a path: a Cmd
    // renders its options before its positionals, and python wants the script
    // ahead of everything
    let cmd = args.only.iter().fold(
        Cmd::new(&args.python)
            .name("plot")
            .sub(script.to_string_lossy())
            .arg("--out", &out)
            .path(&grid),
        |cmd, name| cmd.arg("--only", name),
    );

    let pipeline = PipelineBuilder::new()
        .step(Step::serial([cmd]))
        .stderr_dir(out.join("stderr"))
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

/// That the interpreter is there and can import matplotlib.
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
