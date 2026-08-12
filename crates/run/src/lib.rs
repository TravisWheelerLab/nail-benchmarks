//! Expand a run matrix from a config, execute it against a list of searches,
//! and record what happened.
//!
//! Callers supply the searches and the tool binary location. Nothing here knows
//! about any particular benchmark.

pub mod config;
pub mod exec;
pub mod table;
pub mod tools;

use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use anyhow::{bail, Context};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

/// How often a parallel execution reports progress.
const REPORT_EVERY: Duration = Duration::from_secs(5);

pub use config::{Config, Run};
pub use exec::{Numa, RunsTable};
pub use table::Runs;
pub use tools::{Asset, Bin, Ctx, Search};

/// Finds the repository root from a benchmark crate's manifest directory.
///
/// Call it as `run::repo(env!("CARGO_MANIFEST_DIR"))`. Benchmark crates live
/// at `<repo>/benchmarks/<name>`, so the root is two levels up.
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
    /// How many searches to keep in flight at once.
    ///
    /// One at a time suits searches that are each large enough to use every
    /// core. Many small searches instead want the parallelism here, with
    /// `threads = 1` on each search.
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
pub fn plan(config: &Config, opts: &Options) -> anyhow::Result<Vec<Run>> {
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
/// into the results table. Runs are the outer loop, so a tool's numbers arrive
/// together rather than interleaved shard by shard.
///
/// With `jobs` above one the children compete for the same cores, so the
/// per-row timings become a record of what ran rather than a measurement of it.
pub fn execute(
    config: &Config,
    runs: &[Run],
    searches: &[Search],
    ctx: &Ctx,
    jobs: usize,
) -> anyhow::Result<()> {
    if searches.is_empty() {
        bail!("no searches to run");
    }

    // labels key the prep set and the scratch directories, so a collision
    // would have two searches clearing each other's working directory
    let labels: HashSet<&str> = searches.iter().map(|s| s.label.as_str()).collect();
    if labels.len() != searches.len() {
        bail!("searches must have distinct labels");
    }

    let jobs = jobs.max(1);
    std::fs::create_dir_all(ctx.log_dir())?;

    let table = Mutex::new(RunsTable::create_at(&ctx.runs_table, config.sweep_columns())?);
    let total = runs.len() * searches.len();
    let done = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let errored = Mutex::new(Vec::<String>::new());

    // what: the set of (tool, search) pairs whose prep has been claimed. these
    //       are tracked across the whole execution rather than per outer
    //       iteration.
    //
    // why:  a later run then reuses the databases an earlier one built.
    //
    // note: nothing waits on a claim. a thread that loses one goes straight to
    //       run(), so this holds up only while no two searches in a run share
    //       a label. tools guard their own shared artifacts with build_once.
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
