//! Bringing in result tables that were produced somewhere other than here.
//!
//! Not every search is run by the harness, and that is deliberate: a full
//! hmmsearch over every MGnify shard belongs on a cluster, not on the machine
//! doing the analysis. What comes back is the same tables a pipeline here
//! would have written, timed by whatever `time` that cluster had instead of by
//! `michi`.
//!
//! Nothing is copied. The tables are found where they already are, under
//! `outputs/<pipeline>/results/`, and the only thing written is the
//! `ledger.tbl` beside them -- a result set large enough to be worth running
//! elsewhere is too large to keep a second copy of. A pipeline run here writes
//! its own ledger when it finishes, so this is the one case that needs a
//! command of its own.
//!
//! What a run was is read out of its filename, since that is already the one
//! place every pipeline in this crate records it. See [`Run::parse`] for the
//! grammar.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::Parser;

use util::ledger;
use util::manifest;
use util::time::{self, Timing};

use crate::inputs;

/// The tools a name can identify itself as, which is what says how to read the
/// table.
const TOOLS: [&str; 3] = ["nail", "mmseqs", "hmmer"];

#[derive(Parser, Debug)]
pub struct Args {
    /// A pipeline directory, or the name of one under benchmarks/mgy/outputs/
    #[arg(value_name = "dir|name")]
    pipeline: String,

    /// The tool behind runs whose name doesn't say, as cloud-search's cells don't
    #[arg(long, value_name = "nail|mmseqs|hmmer")]
    tool: Option<String>,

    /// Cover a whole run with one .time file, for a search too long to have
    /// been timed shard by shard. Names a run, or every run without the run=
    #[arg(long, value_name = "[run=]one.time")]
    time: Vec<String>,

    /// Replace an existing ledger.tbl
    #[arg(long)]
    force: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let dir = crate::parse::pipeline(&args.pipeline)?;
    let out = ledger::path(&dir);

    // a pipeline that ran here writes its own ledger, and overwriting it would
    // throw away timings michi measured in favour of ones a shell reported
    ensure!(
        args.force || !out.exists(),
        "{} already exists; pass --force to replace it",
        out.display()
    );

    let ledger = ledger::Ledger::from_rows(rows(
        &dir.join("results"),
        args.tool.as_deref(),
        &args.time,
    )?);
    ledger.write(&out)?;

    println!("wrote {}", out.display());
    for column in ledger.columns()? {
        println!(
            "  {:<24} {:<7} {:>3} shard(s) {:>10.2}s",
            column.name,
            column.tool,
            column.shards.len(),
            column.wall_s
        );
    }

    Ok(())
}

/// Every row the result tables under `results/` and their `.time` files add
/// up to.
///
/// `tool` names the tool behind runs whose own name does not, as
/// cloud-search's cells do not. `time` charges one measurement to a whole run,
/// as `[run=]file`, for a search too long to have been timed shard by shard.
fn rows(results: &Path, tool: Option<&str>, time: &[String]) -> anyhow::Result<Vec<ledger::Row>> {
    ensure!(
        results.is_dir(),
        "no results directory at {}; the tables go there",
        results.display()
    );

    if let Some(tool) = tool {
        ensure!(
            TOOLS.contains(&tool),
            "--tool {tool:?} is not one of {}",
            TOOLS.join(", ")
        );
    }

    let extrapolate = Extrapolate::parse(time)?;

    let mut found = collect(results, tool)?;
    ensure!(
        !found.is_empty(),
        "no result tables in {}; expected <name>.<shard>.tbl",
        results.display()
    );

    for run in &mut found {
        run.check_domains(results)?;
        run.check_targets()?;
        run.read_timings(results, extrapolate.get(&run.name))?;
    }

    Ok(found.iter().flat_map(Run::rows).collect())
}

/// One run's results: what it was, and one entry per shard it covered.
struct Run {
    name: String,
    tool: String,

    /// Whatever the name said past the tool, which is what tells one run of a
    /// tool apart from another.
    params: Vec<(String, String)>,

    shards: Vec<Shard>,

    /// One measurement covering every shard at once, for a search timed as a
    /// whole rather than shard by shard.
    whole: Option<Timing>,
}

struct Shard {
    /// As written, since it travels through the manifest as a string.
    name: String,
    timing: Option<Timing>,
}

impl Run {
    /// What a results filename says about the run that wrote it.
    ///
    /// Every pipeline here names its tables `<run>.<shard>.tbl` through
    /// [`manifest::table_path`], and the run half is built out of `-`-joined
    /// segments. A segment is one of three things:
    ///
    /// ```text
    /// nail    a tool, which settles how the table is read
    /// s12.0   a setting and its value, which becomes a field
    /// rev     anything else, which is part of the name and nothing more
    /// ```
    ///
    /// So `nail-s12.0` is nail at `s=12.0`, `A2.0-B4.0` is `A=2.0` and
    /// `B=4.0` with no tool of its own, and `hmmer` is bare. A name that
    /// identifies no tool needs `--tool`, which is the case for
    /// cloud-search's cells.
    ///
    /// The shard is split off at the last `.` rather than the first, because a
    /// setting's value can hold one: `nail-s12.0.5.tbl` is shard 5 of
    /// `nail-s12.0`, not shard 0 of `nail-s12`.
    fn parse(stem: &str, fallback: Option<&str>) -> anyhow::Result<(Run, String)> {
        let (name, shard) = stem
            .rsplit_once('.')
            .with_context(|| format!("{stem:?} has no shard; expected <name>.<shard>"))?;

        ensure!(
            shard.parse::<usize>().is_ok(),
            "{stem:?} ends in {shard:?}, which is not a shard number"
        );

        let mut tool = None;
        let mut params = Vec::new();

        for segment in name.split('-') {
            if TOOLS.contains(&segment) {
                tool = Some(segment.to_string());
                continue;
            }

            if let Some(param) = setting(segment) {
                params.push(param);
            }
        }

        let tool = tool
            .or_else(|| fallback.map(str::to_string))
            .with_context(|| {
                format!("{name:?} names no tool; pass --tool to say what produced it")
            })?;

        Ok((
            Run {
                name: name.to_string(),
                tool,
                params,
                shards: Vec::new(),
                whole: None,
            },
            shard.to_string(),
        ))
    }

    /// hmmer breaks a hit into domains and `parse` reads that breakdown back,
    /// so a run of it without one is not a run `parse` can use.
    fn check_domains(&self, results: &Path) -> anyhow::Result<()> {
        if self.tool != "hmmer" {
            return Ok(());
        }

        let missing: Vec<PathBuf> = self
            .shards
            .iter()
            .map(|shard| manifest::dom_path(results, &self.name, &shard.name))
            .filter(|path| !path.is_file())
            .collect();

        if missing.is_empty() {
            return Ok(());
        }

        let listed: Vec<String> = missing
            .iter()
            .take(5)
            .map(|p| format!("\n  {}", p.display()))
            .collect();

        let rest = match missing.len() > 5 {
            true => format!("\n  ... {} more", missing.len() - 5),
            false => String::new(),
        };

        bail!(
            "run {:?} is missing {} domain table(s):{}{}\n\n\
             hmmsearch needs --domtblout for this run",
            self.name,
            missing.len(),
            listed.join(""),
            rest
        );
    }

    /// `parse` measures a table against the shard it was searched, so a shard
    /// with nothing to measure against fails there rather than here -- long
    /// after the tables were read.
    fn check_targets(&self) -> anyhow::Result<()> {
        let missing: Vec<&str> = self
            .shards
            .iter()
            .map(|shard| shard.name.as_str())
            .filter(|shard| !inputs::fixed::shard(shard).is_file())
            .collect();

        ensure!(
            missing.is_empty(),
            "run {:?} covers {} shard(s) with no target in {}: {}\n\n\
             `mgy build fixed --shards N` has to have made the same shards",
            self.name,
            missing.len(),
            inputs::fixed::targets().display(),
            missing.join(", ")
        );

        Ok(())
    }

    /// Reads each shard's `.time`, or the one `extrapolate` names.
    ///
    /// Extrapolating is for a search whose whole point is that it was too
    /// expensive to run here: one measurement covers the lot. It is recorded
    /// as one row against every shard rather than copied onto each of them,
    /// so nothing downstream can add it up as though each shard had been
    /// timed.
    fn read_timings(&mut self, results: &Path, extrapolate: Option<&Path>) -> anyhow::Result<()> {
        if let Some(path) = extrapolate {
            self.whole = Some(read_timing(path)?);
            return Ok(());
        }

        for shard in &mut self.shards {
            let path = results.join(format!("{}.{}.time", self.name, shard.name));

            ensure!(
                path.is_file(),
                "no {}\n\n\
                 every table wants a .time beside it, or --time {}=<file> to \
                 cover every shard with one",
                path.display(),
                self.name
            );

            shard.timing = Some(read_timing(&path)?);
        }

        Ok(())
    }

    /// One row per shard, plus the whole-run row where there is one.
    fn rows(&self) -> Vec<ledger::Row> {
        let params: BTreeMap<String, String> = self.params.iter().cloned().collect();
        let row = |shard: &str, wall_s: Option<f64>| ledger::Row {
            name: self.name.clone(),
            tool: self.tool.clone(),
            shard: shard.to_string(),
            stage: String::new(),
            params: params.clone(),
            wall_s,
        };

        let mut out: Vec<ledger::Row> = self
            .shards
            .iter()
            .map(|shard| row(&shard.name, shard.timing.as_ref().map(|t| t.wall_s)))
            .collect();

        if let Some(whole) = &self.whole {
            out.push(row(ledger::EVERY_SHARD, Some(whole.wall_s)));
        }

        out
    }
}

/// One `.time` file, refusing one that recorded a failure.
//
// only GNU `time -v` records a status at all. where there is one, a search
// that failed is not a search to record: its table is half-written, and a
// half-written one reads as a run that simply found less
fn read_timing(path: &Path) -> anyhow::Result<Timing> {
    let timing = time::read(path)?;

    if let Some(exit) = timing.exit.filter(|exit| *exit != 0) {
        bail!(
            "{} recorded exit status {exit}; that search did not finish",
            path.display()
        );
    }

    Ok(timing)
}

/// One `-`-separated segment read as a setting: a name, then a value that has
/// to look like a number.
///
/// `s12.0` is a setting and `rev` is not, which is what keeps a direction or a
/// mode out of the parameter columns.
fn setting(segment: &str) -> Option<(String, String)> {
    let at = segment.find(|c: char| !c.is_ascii_alphabetic())?;
    let (key, value) = segment.split_at(at);

    if key.is_empty() || value.parse::<f64>().is_err() {
        return None;
    }

    Some((key.to_string(), value.to_string()))
}

/// Every run in `results/`, grouped by name.
///
/// Ordered by tool and then by name, which is the order the pipelines here
/// declare their runs in -- so an installed manifest puts the columns of a
/// scores table where a native one would, rather than wherever the directory
/// happened to list them.
fn collect(results: &Path, fallback: Option<&str>) -> anyhow::Result<Vec<Run>> {
    let mut runs: BTreeMap<String, Run> = BTreeMap::new();

    let dir = std::fs::read_dir(results)
        .with_context(|| format!("failed to read {}", results.display()))?;

    for entry in dir {
        let path = entry?.path();

        // the domain tables and the seed lists live here too, and neither is a
        // run of its own
        if path.extension().is_none_or(|ext| ext != "tbl") {
            continue;
        }

        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .with_context(|| format!("{} has no readable name", path.display()))?;

        let (run, shard) = Run::parse(stem, fallback)?;

        runs.entry(run.name.clone())
            .or_insert(run)
            .shards
            .push(Shard {
                name: shard,
                timing: None,
            });
    }

    // numeric, so shard 10 does not sort between 1 and 2
    for run in runs.values_mut() {
        run.shards.sort_by_key(|shard| {
            shard
                .name
                .parse::<usize>()
                .expect("parse checked the shard is a number")
        });
    }

    let mut runs: Vec<Run> = runs.into_values().collect();
    runs.sort_by(|x, y| {
        tool_rank(x)
            .cmp(&tool_rank(y))
            .then_with(|| swept(x).partial_cmp(&swept(y)).unwrap_or(Ordering::Equal))
            .then_with(|| x.name.cmp(&y.name))
    });

    Ok(runs)
}

fn tool_rank(run: &Run) -> Option<usize> {
    TOOLS.iter().position(|tool| *tool == run.tool)
}

/// What a run swept, for ordering one against another.
///
/// The values compare as numbers rather than as text, so a sweep comes out in
/// the order it was run rather than with `A10.0` in front of `A2.0`. Sorted by
/// name first, so two runs are compared knob by knob whatever order their
/// names listed them in.
fn swept(run: &Run) -> Vec<(&str, f64)> {
    let mut out: Vec<(&str, f64)> = run
        .params
        .iter()
        .map(|(key, value)| {
            // setting() only makes a param out of a value that parses
            (key.as_str(), value.parse().unwrap_or(f64::NAN))
        })
        .collect();

    out.sort_by(|x, y| x.0.cmp(y.0));
    out
}

/// Which `.time` file stands in for a run's own, from the `--time` arguments.
#[derive(Default)]
struct Extrapolate {
    /// For every run that was not named, which is the whole pipeline when one
    /// machine timed one shard of it.
    every: Option<PathBuf>,
    by_run: BTreeMap<String, PathBuf>,
}

impl Extrapolate {
    fn parse(args: &[String]) -> anyhow::Result<Extrapolate> {
        let mut out = Extrapolate::default();

        for arg in args {
            let (run, path) = match arg.split_once('=') {
                Some((run, path)) => (Some(run), path),
                None => (None, arg.as_str()),
            };

            let path = PathBuf::from(path);
            ensure!(path.is_file(), "no {} for --time", path.display());

            match run {
                Some(run) => {
                    out.by_run.insert(run.to_string(), path);
                }
                None => {
                    ensure!(
                        out.every.is_none(),
                        "--time was given twice with no run to charge it to"
                    );
                    out.every = Some(path);
                }
            }
        }

        Ok(out)
    }

    fn get(&self, run: &str) -> Option<&Path> {
        self.by_run
            .get(run)
            .map(PathBuf::as_path)
            .or(self.every.as_deref())
    }
}
