//! Draws the figures, by handing what `parse` wrote to matplotlib.
//!
//! The drawing is python because matplotlib is what the other benchmarks plot
//! with. It goes through the pipeline like everything else, which is what gets
//! it a `--dry-run`, its stderr kept on failure, and a line in the progress
//! output.
//!
//! Each script takes its input and its output pdf and nothing else, so a
//! figure here is a name, a script, and which of `parse`'s files it reads.

use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, bail};
use clap::Parser;

use pail::{Cmd, PipelineBuilder, Progress, Step};

use crate::inputs::{self, Inputs};

/// A figure: which script draws it, and which of `parse`'s files it reads.
///
/// `plot_params` and `plot_threads` are not here. They read a points file and
/// a directory of thread sweeps that nothing in this benchmark produces --
/// there is no `--threads` axis in the run, and no producer for the points
/// survived. Wiring them up would mean deciding what they should read.
struct Figure {
    name: &'static str,
    script: &'static str,
    inputs: &'static [&'static str],
}

const FIGURES: &[Figure] = &[
    Figure {
        name: "roc",
        script: "plot_roc.py",
        inputs: &["roc.txt"],
    },
    Figure {
        name: "pid",
        script: "plot_pid.py",
        inputs: &["pid.txt"],
    },
    Figure {
        name: "time",
        script: "plot_time.py",
        inputs: &["time.txt"],
    },
    Figure {
        name: "cells",
        script: "plot_cell_frac.py",
        inputs: &["cells.true.txt", "cells.decoy.txt"],
    },
    Figure {
        name: "score",
        script: "plot_score.py",
        inputs: &["score.txt"],
    },
];

#[derive(Parser, Debug)]
pub struct Args {
    /// Which input set's figures to draw, naming runs/<size>/.
    #[arg(short, long, default_value = "toy")]
    size: String,

    /// Where `parse` wrote its tables, and where the pdfs go. Defaults to
    /// figures/ beside the run.
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,

    /// Draw only this figure. Repeatable; every one it has the input for by
    /// default.
    #[arg(long, value_name = "NAME")]
    only: Vec<String>,

    /// The interpreter to run the scripts with.
    #[arg(long, default_value = "python3", value_name = "python")]
    python: String,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let names: Vec<&'static str> = FIGURES.iter().map(|f| f.name).collect();
    for only in &args.only {
        if !names.contains(&only.as_str()) {
            bail!("no figure named {only:?}; there is {}", names.join(", "));
        }
    }

    let dir = args
        .out
        .unwrap_or_else(|| Inputs::new(&args.size).run_dir().join("figures"));

    if !dir.is_dir() {
        bail!(
            "no {}; run `pct-id parse recall --size {}` first",
            dir.display(),
            args.size
        );
    }

    let scripts = inputs::dir().join("scripts");

    // checked here rather than left to the pipeline, because a missing
    // matplotlib comes back as a python traceback in a stderr file rather than
    // as anything that reads like an answer
    matplotlib(&args.python)?;

    let mut pl = PipelineBuilder::new();
    let mut drawn = 0usize;

    for figure in FIGURES {
        if !args.only.is_empty() && !args.only.iter().any(|o| o == figure.name) {
            continue;
        }

        // a figure whose analysis has not been run is skipped rather than
        // failed: `parse recall` and `parse cells` write different files, and
        // asking for everything should not mean running everything first
        let missing: Vec<&str> = figure
            .inputs
            .iter()
            .copied()
            .filter(|f| !dir.join(f).is_file())
            .collect();

        if !missing.is_empty() {
            eprintln!(
                "skipping {}: no {} in {}",
                figure.name,
                missing.join(", "),
                dir.display()
            );
            continue;
        }

        let script = scripts.join(figure.script);
        if !script.is_file() {
            bail!("no plot script at {}", script.display());
        }

        // the script goes in the subcommand slot rather than in a path: a Cmd
        // renders its options before its positionals, and python wants the
        // script ahead of everything
        let cmd = figure
            .inputs
            .iter()
            .fold(
                Cmd::new(&args.python)
                    .name(figure.name)
                    .sub(script.to_string_lossy()),
                |cmd, input| cmd.path(dir.join(input)),
            )
            .path(dir.join(format!("{}.pdf", figure.name)));

        pl = pl.step(Step::serial([cmd]).name(figure.name));
        drawn += 1;
    }

    if drawn == 0 {
        bail!(
            "nothing to draw: no analysis output in {}. run `pct-id parse recall` first",
            dir.display()
        );
    }

    let pipeline = pl
        .stderr_dir(dir.join("stderr"))
        .sink(Progress::new())
        .build()
        .context("failed to build the figures")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()?;

    println!("\nfigures in {}", dir.display());
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
