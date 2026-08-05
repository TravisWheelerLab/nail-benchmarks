//! The generic half of the benchmark harness: expand a run matrix from a
//! config, execute it against a list of searches, and record what happened.
//!
//! This library knows nothing about any particular benchmark. Callers supply
//! the searches ([`Search`]) and the tool binary location; how those were
//! constructed is the benchmark's business.

pub mod config;
pub mod exec;
pub mod split;
pub mod tools;

use std::collections::HashSet;

use anyhow::{bail, Context, Result};

pub use config::{Config, Run};
pub use exec::{Numa, RunsTable};
pub use tools::{Asset, Bin, Ctx, Search};

/// Options that override what the config declares.
#[derive(Default, Debug)]
pub struct Options {
    pub filter: Option<String>,
    pub threads: Option<usize>,
    pub numa_node: Option<usize>,
    pub dry_run: bool,
}

/// Expand a config's run matrix, applying overrides and any name filter.
pub fn plan(config: &Config, opts: &Options) -> Result<Vec<Run>> {
    let mut runs = config.expand()?;

    if let Some(pattern) = &opts.filter {
        let pattern = glob::Pattern::new(pattern)
            .with_context(|| format!("invalid filter glob {pattern:?}"))?;
        runs.retain(|r| pattern.matches(&r.name));

        if runs.is_empty() {
            bail!("filter {pattern:?} matched no runs");
        }
    }

    if let Some(threads) = opts.threads {
        for run in &mut runs {
            run.threads = threads;
        }
    }

    Ok(runs)
}

/// Execute every run against every search, writing one row per (run, search)
/// into the results table.
///
/// Searches are the outer loop so that per-target setup — an mmseqs target db,
/// a blast db — is built once and then reused by every run against it.
pub fn execute(
    config: &Config,
    runs: &[Run],
    searches: &[Search],
    ctx: &Ctx,
) -> Result<()> {
    if searches.is_empty() {
        bail!("no searches to run");
    }

    let mut table = RunsTable::create(&ctx.results, config.sweep_columns())?;
    let total = runs.len() * searches.len();
    let mut done = 0usize;
    let mut failed = 0usize;

    for search in searches {
        // (tool, search) pairs already prepared; tools additionally skip work
        // whose output exists, so a shared query set converts only once
        let mut prepped: HashSet<String> = HashSet::new();

        for run in runs {
            let tool = tools::get(&run.tool)?;

            if prepped.insert(run.tool.clone()) {
                tool.prep(ctx, search)
                    .with_context(|| format!("prep failed for tool {:?}", run.tool))?;
                std::fs::remove_file(ctx.log_dir().join(format!("prep-{}.err", run.tool))).ok();
            }

            done += 1;
            let label = if search.label.is_empty() {
                run.name.clone()
            } else {
                format!("{} [{}]", run.name, search.label)
            };
            println!("[{done}/{total}] {label} | {}", run.args.join(" "));

            let outcome = tool
                .run(ctx, search, run)
                .with_context(|| format!("run {:?} failed to execute", run.name))?;

            let log = ctx.log_path(run, search);
            if outcome.timing.exit != 0 {
                failed += 1;
                eprintln!(
                    "  FAILED: exit {} after {:.2}s",
                    outcome.timing.exit, outcome.timing.wall_s
                );
                for line in tail(&log, 10) {
                    eprintln!("  | {line}");
                }
                eprintln!("  full stderr: {}", log.display());
            } else {
                println!(
                    "  {:.2}s wall, {:.2}s cpu, {} kb peak",
                    outcome.timing.wall_s,
                    outcome.timing.user_s + outcome.timing.sys_s,
                    outcome.timing.max_rss_kb
                );
                // tools are chatty on success; the log is only worth keeping
                // when something went wrong
                std::fs::remove_file(&log).ok();
            }

            table.append(run, &search.display(), &outcome.timing, &outcome.cmd)?;
        }
    }

    // leftover logs belong to failed runs; drop the directory when it is empty
    std::fs::remove_dir(ctx.log_dir()).ok();

    println!("\nwrote {}", table.path().display());

    if failed > 0 {
        bail!("{failed} of {total} runs exited nonzero");
    }

    Ok(())
}

/// Print the expanded matrix without executing it.
pub fn describe(runs: &[Run], searches: &[Search]) {
    for run in runs {
        println!("{:<34} {:<8} {}", run.name, run.tool, run.args.join(" "));
    }
    println!(
        "\n{} runs x {} searches = {} executions",
        runs.len(),
        searches.len(),
        runs.len() * searches.len()
    );
}

/// Last `n` non-empty lines of a file, for surfacing why a run failed.
fn tail(path: &std::path::Path, n: usize) -> Vec<String> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };

    let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
    lines
        .iter()
        .skip(lines.len().saturating_sub(n))
        .map(|l| l.to_string())
        .collect()
}
