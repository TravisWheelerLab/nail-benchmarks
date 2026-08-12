//! Expand a run matrix from a config, execute it against a list of searches,
//! and record what happened.
//!
//! Callers supply the searches and the tool binary location. Nothing here knows
//! about any particular benchmark.
//!
//! There are two executors. [`measure`] runs one search at a time, so a row's
//! timing is a measurement of that search alone. [`batch`] runs several at once,
//! for stages where the output tables are the point and the clock is not.

pub mod config;
pub mod exec;
pub mod table;
pub mod tools;

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use anyhow::{bail, Context};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

pub use config::{Config, Run};
pub use exec::{Numa, Timing};
pub use table::Runs;
pub use tools::Tool;

use exec::Cmd;
use tools::{After, Shape};

/// How often a batch execution reports progress.
const REPORT_EVERY: Duration = Duration::from_secs(5);

/// Finds the repository root from a benchmark crate's manifest directory.
///
/// Call it as `run::repo(env!("CARGO_MANIFEST_DIR"))`. Benchmark crates live
/// at `<repo>/benchmarks/<name>`, so the root is two levels up.
pub fn repo(manifest_dir: &str) -> PathBuf {
    Path::new(manifest_dir)
        .parent()
        .and_then(Path::parent)
        .expect("benchmark crates live at <repo>/benchmarks/<name>")
        .to_path_buf()
}

// ------------------------------------------------------------------ searches

/// One unit of work: a set of query artifacts searched against one target.
///
/// The benchmark builds this list, so the runner never has to guess at a
/// directory layout. Which artifacts a tool needs depends on its query mode,
/// and a tool asked for one that is absent fails with a message naming it.
#[derive(Clone, Debug)]
pub struct Search {
    /// Distinguishes this search's outputs from the others in the same run.
    /// Empty when a benchmark has only one search.
    pub label: String,
    pub target: PathBuf,
    /// HMMER3 profiles.
    pub hmm: Option<PathBuf>,
    /// Stockholm alignments.
    pub sto: Option<PathBuf>,
    /// Unaligned query sequences.
    pub fasta: Option<PathBuf>,
    /// Directory of per-family aligned fasta.
    pub afa: Option<PathBuf>,
}

impl Search {
    pub fn new(label: impl Into<String>, target: impl Into<PathBuf>) -> Self {
        Search {
            label: label.into(),
            target: target.into(),
            hmm: None,
            sto: None,
            fasta: None,
            afa: None,
        }
    }

    pub fn with_hmm(mut self, path: impl Into<PathBuf>) -> Self {
        self.hmm = Some(path.into());
        self
    }

    pub fn with_sto(mut self, path: impl Into<PathBuf>) -> Self {
        self.sto = Some(path.into());
        self
    }

    pub fn with_fasta(mut self, path: impl Into<PathBuf>) -> Self {
        self.fasta = Some(path.into());
        self
    }

    pub fn with_afa(mut self, path: impl Into<PathBuf>) -> Self {
        self.afa = Some(path.into());
        self
    }
}

// --------------------------------------------------------------------- paths

/// Where a benchmark keeps its binaries, scratch space and results.
pub struct Paths {
    /// Directory holding the tool binaries, usually `<repo>/tools/bin`.
    pub bin: PathBuf,
    pub tmp: PathBuf,
    pub results: PathBuf,
    /// Where the runs table is written. Usually inside `results`, but a
    /// benchmark may keep it above so clearing results does not remove it.
    pub runs_table: PathBuf,
    pub numa: Option<Numa>,
}

/// Locate a tool binary, naming the make target that produces it when it is
/// not there.
///
/// Canonicalized, so the `cmd` column of the runs table shows a path you can
/// paste into a shell.
pub fn tool(bin_dir: impl AsRef<Path>, name: &str) -> anyhow::Result<PathBuf> {
    let path = bin_dir.as_ref().join(name);
    if !path.exists() {
        bail!(
            "missing tool binary {}; run `make {}` from the repo root",
            path.display(),
            make_target(name)
        );
    }
    Ok(path.canonicalize().unwrap_or(path))
}

impl Paths {
    pub fn tool(&self, name: &str) -> anyhow::Result<PathBuf> {
        tool(&self.bin, name)
    }

    pub fn log_dir(&self) -> PathBuf {
        self.results.join(".logs")
    }

    /// Output path for one run of one search.
    pub fn out(&self, run: &Run, search: &Search, ext: &str) -> PathBuf {
        self.results.join(format!("{}.{ext}", stem(run, search)))
    }

    pub fn log(&self, run: &Run, search: &Search) -> PathBuf {
        self.log_dir().join(format!("{}.err", stem(run, search)))
    }

    /// Working directory for one (run, search), emptied before the search runs.
    pub fn scratch(&self, tool: Tool, run: &Run, search: &Search) -> PathBuf {
        self.tmp.join(tool.name()).join(stem(run, search))
    }
}

/// What names one (run, search) pair's files: outputs, logs and scratch all
/// share it, so a row in the runs table and the files it refers to agree.
fn stem(run: &Run, search: &Search) -> String {
    if search.label.is_empty() {
        run.name.clone()
    } else {
        format!("{}.{}", run.name, search.label)
    }
}

fn make_target(bin: &str) -> &'static str {
    match bin {
        "nail" => "nail",
        "hmmsearch" | "phmmer" | "hmmbuild" | "hmmemit" | "esl-seqstat" | "create-profmark" => {
            "hmmer"
        }
        "mmseqs" => "mmseqs",
        "blastp" | "psiblast" | "makeblastdb" => "blast",
        "lastal" | "lastdb" => "last",
        "diamond" => "diamond",
        _ => "setup",
    }
}

// -------------------------------------------------------------------- config

/// Options that override what the config declares.
#[derive(Debug, Default)]
pub struct Options {
    pub filter: Option<String>,
    pub threads: Option<usize>,
    pub numa_node: Option<usize>,
    pub dry_run: bool,
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

// ----------------------------------------------------------------- executors

/// Execute every run against every search, one at a time, writing one row per
/// pair into the results table.
///
/// Runs are the outer loop, so a tool's numbers arrive together rather than
/// interleaved shard by shard. Nothing else is competing for the machine, so
/// each row's timing measures that search alone.
pub fn measure(
    config: &Config,
    runs: &[Run],
    searches: &[Search],
    paths: &Paths,
) -> anyhow::Result<()> {
    let matrix = prepare(runs, searches, paths)?;
    let mut table = table::Writer::create(&paths.runs_table, config.sweep_columns())?;
    let report = Report::new(runs.len() * searches.len());

    for (tool, run) in &matrix {
        for search in searches {
            let label = describe_one(run, search);
            println!("[{}/{}] {label} | {}", report.next(), report.total, run.args.join(" "));

            match execute_one(*tool, run, search, paths) {
                Err(e) => report.error(format!("{label}: {e:#}")),
                Ok((timing, cmd)) => {
                    if timing.exit == 0 {
                        println!(
                            "  {:.2}s wall, {:.2}s cpu, {} kb peak",
                            timing.wall_s,
                            timing.user_s + timing.sys_s,
                            timing.max_rss_kb
                        );
                    } else {
                        report.failure(&label, &timing, &paths.log(run, search));
                    }
                    record(&mut table, run, search, &timing, &cmd, &label);
                }
            }
        }
    }

    finish(table, paths, report)
}

/// Execute every run against every search, `jobs` searches at a time.
///
/// The children compete for the same cores, so the per-row timings become a
/// record of what ran rather than a measurement of it. Use this where the hit
/// tables are what matters and the searches are each too small to fill the
/// machine on their own.
pub fn batch(
    config: &Config,
    runs: &[Run],
    searches: &[Search],
    paths: &Paths,
    jobs: usize,
) -> anyhow::Result<()> {
    let jobs = jobs.max(1);
    let matrix = prepare(runs, searches, paths)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .context("failed to build the executor thread pool")?;

    let table = Mutex::new(table::Writer::create(
        &paths.runs_table,
        config.sweep_columns(),
    )?);
    let report = Report::new(runs.len() * searches.len());
    let last_report = Mutex::new(Instant::now());

    for (tool, run) in &matrix {
        pool.install(|| {
            searches.par_iter().for_each(|search| {
                let label = describe_one(run, search);
                let n = report.next();

                match execute_one(*tool, run, search, paths) {
                    Err(e) => report.error(format!("{label}: {e:#}")),
                    Ok((timing, cmd)) => {
                        if timing.exit != 0 {
                            report.failure(&label, &timing, &paths.log(run, search));
                        }
                        let mut table = table.lock().expect("table mutex poisoned");
                        record(&mut table, run, search, &timing, &cmd, &label);
                    }
                }

                // per-execution chatter is noise at this width; report on a
                // timer instead
                let mut last = last_report.lock().expect("report mutex poisoned");
                if last.elapsed() >= REPORT_EVERY {
                    *last = Instant::now();
                    report.progress(n, &run.name);
                }
            })
        });
    }

    finish(
        table.into_inner().expect("table mutex poisoned"),
        paths,
        report,
    )
}

// ------------------------------------------------------------------ internals

/// Check the matrix, resolve its tools, and build every database it needs.
///
/// Both executors run this before measuring anything, so a bad config or a
/// missing binary surfaces before any time is spent.
fn prepare(
    runs: &[Run],
    searches: &[Search],
    paths: &Paths,
) -> anyhow::Result<Vec<(Tool, Run)>> {
    if searches.is_empty() {
        bail!("no searches to run");
    }

    // labels name the outputs, the logs and the scratch directories, so a
    // collision would have two searches overwriting each other everywhere
    let labels: HashSet<&str> = searches.iter().map(|s| s.label.as_str()).collect();
    if labels.len() != searches.len() {
        bail!("searches must have distinct labels");
    }

    let runs: Vec<(Tool, Run)> = runs
        .iter()
        .map(|r| Ok((Tool::parse(&r.tool)?, r.clone())))
        .collect::<anyhow::Result<_>>()?;

    std::fs::create_dir_all(paths.log_dir())?;
    build_databases(&runs, searches, paths)?;

    Ok(runs)
}

/// Build every database the matrix needs, before anything is measured.
///
/// This runs to completion on one thread, so two searches deriving the same
/// database path cannot race to build it and leave one that is neither.
fn build_databases(
    runs: &[(Tool, Run)],
    searches: &[Search],
    paths: &Paths,
) -> anyhow::Result<()> {
    let mut tools: Vec<Tool> = Vec::new();
    for (tool, _) in runs {
        if !tools.contains(tool) {
            tools.push(*tool);
        }
    }

    let mut done: HashSet<PathBuf> = HashSet::new();

    for tool in tools {
        let log = paths.log_dir().join(format!("prep-{}.err", tool.name()));
        let mut built = 0;

        for search in searches {
            for step in tool.setup(search, paths)? {
                if !done.insert(step.marker.clone()) || step.marker.exists() {
                    continue;
                }

                for dir in &step.dirs {
                    std::fs::create_dir_all(dir)
                        .with_context(|| format!("failed to create {}", dir.display()))?;
                }
                for cmd in &step.cmds {
                    let cmd = cmd.clone().stderr_to(&log);
                    exec::check(&cmd, paths.numa.as_ref(), &step.what)
                        .with_context(|| format!("setup failed for tool {}", tool.name()))?;
                }
                built += 1;
            }
        }

        if built > 0 {
            println!("built {built} {} databases", tool.name());
            std::fs::remove_file(&log).ok();
        }
    }

    Ok(())
}

/// Run one (run, search) pair and hand back its timing and command line.
fn execute_one(
    tool: Tool,
    run: &Run,
    search: &Search,
    paths: &Paths,
) -> anyhow::Result<(Timing, String)> {
    let numa = paths.numa.as_ref();

    if tool.uses_scratch() {
        let dir = paths.scratch(tool, run, search);
        if dir.exists() {
            std::fs::remove_dir_all(&dir)?;
        }
        std::fs::create_dir_all(&dir)?;
    }

    let work = tool
        .work(run, search, paths)
        .with_context(|| format!("could not plan run {:?}", run.name))?;

    for dir in &work.dirs {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("failed to create {}", dir.display()))?;
    }

    let (timing, cmd) = match &work.search {
        Shape::One(cmd) => (exec::run(cmd, numa)?, exec::render(cmd, numa)),
        Shape::Together(cmds) => (exec::run_together(cmds, numa)?, repeated(cmds, numa)),
        Shape::Each(cmds) => (exec::run_each(cmds, numa)?, repeated(cmds, numa)),
    };

    // a failed search leaves nothing worth converting or concatenating, and its
    // exit code is already on its way to the table
    if timing.exit == 0 {
        for step in &work.after {
            apply(step, numa)?;
        }
        // tools are chatty on success; the log is only worth keeping when
        // something went wrong
        std::fs::remove_file(paths.log(run, search)).ok();
    }

    Ok((timing, cmd))
}

fn repeated(cmds: &[Cmd], numa: Option<&Numa>) -> String {
    format!("[{} x] {}", cmds.len(), exec::render(&cmds[0], numa))
}

fn apply(step: &After, numa: Option<&Numa>) -> anyhow::Result<()> {
    match step {
        After::Concat { parts, into } => {
            if let Some(dir) = into.parent() {
                std::fs::create_dir_all(dir)?;
            }
            let mut dst = std::fs::File::create(into)
                .with_context(|| format!("failed to create {}", into.display()))?;

            for part in parts {
                let mut src = std::fs::File::open(part)
                    .with_context(|| format!("failed to open {}", part.display()))?;
                std::io::copy(&mut src, &mut dst)?;
            }
            Ok(())
        }
        After::Move { from, to } => {
            if from.exists() {
                std::fs::rename(from, to)?;
            }
            Ok(())
        }
        After::Run { cmd, what } => exec::check(cmd, numa, what),
        After::Remove(path) => {
            std::fs::remove_dir_all(path).ok();
            Ok(())
        }
    }
}

fn describe_one(run: &Run, search: &Search) -> String {
    if search.label.is_empty() {
        run.name.clone()
    } else {
        format!("{} [{}]", run.name, search.label)
    }
}

fn record(
    table: &mut table::Writer,
    run: &Run,
    search: &Search,
    timing: &Timing,
    cmd: &str,
    label: &str,
) {
    if let Err(e) = table.append(run, &search.label, timing, cmd) {
        eprintln!("  WARNING: could not record {label}: {e:#}");
    }
}

fn finish(mut table: table::Writer, paths: &Paths, report: Report) -> anyhow::Result<()> {
    table.flush()?;

    // leftover logs belong to failed runs; drop the directory when it is empty
    std::fs::remove_dir(paths.log_dir()).ok();

    println!("\nwrote {}", table.path().display());
    report.into_result()
}

/// Progress and failure accounting, shared by the two executors.
struct Report {
    total: usize,
    done: AtomicUsize,
    failed: AtomicUsize,
    errors: Mutex<Vec<String>>,
}

impl Report {
    fn new(total: usize) -> Self {
        Report {
            total,
            done: AtomicUsize::new(0),
            failed: AtomicUsize::new(0),
            errors: Mutex::new(Vec::new()),
        }
    }

    fn next(&self) -> usize {
        self.done.fetch_add(1, Ordering::Relaxed) + 1
    }

    fn error(&self, message: String) {
        self.failed.fetch_add(1, Ordering::Relaxed);
        eprintln!("  ERROR {message}");
        self.errors
            .lock()
            .expect("errors mutex poisoned")
            .push(message);
    }

    fn failure(&self, label: &str, timing: &Timing, log: &Path) {
        self.failed.fetch_add(1, Ordering::Relaxed);
        eprintln!(
            "  FAILED {label}: exit {} after {:.2}s",
            timing.exit, timing.wall_s
        );
        for line in tail(log, 10) {
            eprintln!("  | {line}");
        }
        eprintln!("  full stderr: {}", log.display());
    }

    fn progress(&self, n: usize, name: &str) {
        let bad = self.failed.load(Ordering::Relaxed);
        let note = if bad > 0 {
            format!(", {bad} failed")
        } else {
            String::new()
        };
        println!("[{n}/{}] {name}{note}", self.total);
    }

    fn into_result(self) -> anyhow::Result<()> {
        let failed = self.failed.load(Ordering::Relaxed);
        if failed == 0 {
            return Ok(());
        }

        let errors = self.errors.into_inner().expect("errors mutex poisoned");
        for e in errors.iter().take(5) {
            eprintln!("  {e}");
        }
        if errors.len() > 5 {
            eprintln!("  ... and {} more", errors.len() - 5);
        }

        bail!("{failed} of {} runs failed", self.total)
    }
}

/// Last `n` non-empty lines of a file, for surfacing why a run failed.
fn tail(path: &Path, n: usize) -> Vec<String> {
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
