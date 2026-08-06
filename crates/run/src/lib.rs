//! The generic half of the benchmark harness: expand a run matrix from a
//! config, execute it against a list of searches, and record what happened.
//!
//! This library knows nothing about any particular benchmark. Callers supply
//! the searches ([`Search`]) and the tool binary location; how those were
//! constructed is the benchmark's business.

pub mod config;
pub mod exec;
pub mod table;
pub mod tools;

use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

/// How often a parallel execution reports progress.
const REPORT_EVERY: Duration = Duration::from_secs(5);

pub use config::{Config, Run};
pub use exec::{Numa, RunsTable};
pub use table::Runs;
pub use tools::{Asset, Bin, Ctx, Search};

/// The repository root, derived from a benchmark crate's manifest directory.
///
/// Call as `run::repo(env!("CARGO_MANIFEST_DIR"))`: benchmark crates live at
/// `<repo>/benchmarks/<name>`, so the root is two levels up. This is fixed at
/// compile time, which is fine because the binaries are built and run out of
/// the same checkout.
pub fn repo(manifest_dir: &str) -> std::path::PathBuf {
    std::path::Path::new(manifest_dir)
        .parent()
        .and_then(std::path::Path::parent)
        .expect("benchmark crates live at <repo>/benchmarks/<name>")
        .to_path_buf()
}

/// Options that override what the config declares.
#[derive(Debug)]
pub struct Options {
    pub filter: Option<String>,
    pub threads: Option<usize>,
    pub numa_node: Option<usize>,
    /// How many searches to keep in flight at once. One search at a time is
    /// right when each is large enough to use every core; a calibration pass
    /// made of many small searches wants the parallelism on the outside
    /// instead, with `threads = 1` per search.
    pub jobs: usize,
    pub dry_run: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            filter: None,
            threads: None,
            numa_node: None,
            jobs: 1,
            dry_run: false,
        }
    }
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
/// Runs are the outer loop: one configuration sweeps every search before the
/// next begins, so a tool's numbers arrive together rather than interleaved
/// shard by shard.
///
/// `jobs` is how many searches run concurrently, each with its own `threads`.
/// Above one, the per-row timings stop being benchmark measurements — the
/// children are competing for the same cores — so treat them as a record of
/// what ran rather than as numbers to plot.
pub fn execute(
    config: &Config,
    runs: &[Run],
    searches: &[Search],
    ctx: &Ctx,
    jobs: usize,
) -> Result<()> {
    if searches.is_empty() {
        bail!("no searches to run");
    }

    let jobs = jobs.max(1);
    std::fs::create_dir_all(ctx.log_dir())?;

    let table = Mutex::new(RunsTable::create_at(&ctx.runs_table, config.sweep_columns())?);
    let total = runs.len() * searches.len();
    let done = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let errored = Mutex::new(Vec::<String>::new());

    // (tool, search) pairs already prepared. Tracked across the whole execution
    // rather than per outer iteration, so a later run reuses the databases an
    // earlier one built; tools additionally skip work whose output exists.
    //
    // Claiming a key and doing the prep are not atomic together, which is safe
    // only because runs are sequential and each search appears once per run, so
    // no two threads ever contend for the same key.
    let prepped = Mutex::new(HashSet::<(String, String)>::new());

    // when several searches are in flight, per-execution chatter is noise;
    // report on a timer instead
    let last_report = Mutex::new(Instant::now());

    let pool = (jobs > 1)
        .then(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(jobs)
                .build()
                .context("failed to build the executor thread pool")
        })
        .transpose()?;

    for run in runs {
        let tool = tools::get(&run.tool)?;

        let one = |search: &Search| {
            let claimed = {
                let mut set = prepped.lock().expect("prepped mutex poisoned");
                set.insert((run.tool.clone(), search.label.clone()))
            };

            if claimed {
                if let Err(e) = tool
                    .prep(ctx, search)
                    .with_context(|| format!("prep failed for tool {:?}", run.tool))
                {
                    failed.fetch_add(1, Ordering::Relaxed);
                    errored
                        .lock()
                        .expect("errored mutex poisoned")
                        .push(format!("{e:#}"));
                    return;
                }
                std::fs::remove_file(ctx.log_dir().join(format!("prep-{}.err", run.tool))).ok();
            }

            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            let label = if search.label.is_empty() {
                run.name.clone()
            } else {
                format!("{} [{}]", run.name, search.label)
            };

            if jobs == 1 {
                println!("[{n}/{total}] {label} | {}", run.args.join(" "));
            }

            let outcome = match tool
                .run(ctx, search, run)
                .with_context(|| format!("run {:?} failed to execute", run.name))
            {
                Ok(o) => o,
                Err(e) => {
                    failed.fetch_add(1, Ordering::Relaxed);
                    eprintln!("  ERROR {label}: {e:#}");
                    errored
                        .lock()
                        .expect("errored mutex poisoned")
                        .push(format!("{e:#}"));
                    return;
                }
            };

            let log = ctx.log_path(run, search);
            if outcome.timing.exit != 0 {
                failed.fetch_add(1, Ordering::Relaxed);
                eprintln!(
                    "  FAILED {label}: exit {} after {:.2}s",
                    outcome.timing.exit, outcome.timing.wall_s
                );
                for line in tail(&log, 10) {
                    eprintln!("  | {line}");
                }
                eprintln!("  full stderr: {}", log.display());
            } else {
                if jobs == 1 {
                    println!(
                        "  {:.2}s wall, {:.2}s cpu, {} kb peak",
                        outcome.timing.wall_s,
                        outcome.timing.user_s + outcome.timing.sys_s,
                        outcome.timing.max_rss_kb
                    );
                }
                // tools are chatty on success; the log is only worth keeping
                // when something went wrong
                std::fs::remove_file(&log).ok();
            }

            if let Err(e) = table.lock().expect("table mutex poisoned").append(
                run,
                &search.display(),
                &outcome.timing,
                &outcome.cmd,
            ) {
                eprintln!("  WARNING: could not record {label}: {e:#}");
            }

            if jobs > 1 {
                let mut last = last_report.lock().expect("report mutex poisoned");
                if last.elapsed() >= REPORT_EVERY {
                    *last = Instant::now();
                    let bad = failed.load(Ordering::Relaxed);
                    let note = if bad > 0 {
                        format!(", {bad} failed")
                    } else {
                        String::new()
                    };
                    println!("[{n}/{total}] {}{note}", run.name);
                }
            }
        };

        match &pool {
            Some(pool) => pool.install(|| searches.par_iter().for_each(one)),
            None => searches.iter().for_each(one),
        }
    }

    let mut table = table.into_inner().expect("table mutex poisoned");
    table.flush()?;

    // leftover logs belong to failed runs; drop the directory when it is empty
    std::fs::remove_dir(ctx.log_dir()).ok();

    println!("\nwrote {}", table.path().display());

    let failed = failed.load(Ordering::Relaxed);
    if failed > 0 {
        let errors = errored.into_inner().expect("errored mutex poisoned");
        for e in errors.iter().take(5) {
            eprintln!("  {e}");
        }
        if errors.len() > 5 {
            eprintln!("  ... and {} more", errors.len() - 5);
        }
        bail!("{failed} of {total} runs failed");
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
