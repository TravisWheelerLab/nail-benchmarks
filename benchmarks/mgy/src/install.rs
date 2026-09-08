//! Turning a directory of externally-run result tables into a pipeline.
//!
//! Not every search is run by the harness, and that is deliberate: a full
//! hmmsearch over every MGnify shard belongs on a cluster, not on the machine
//! doing the analysis. What comes back is the same tables a pipeline here
//! would have written, timed by whatever `time` that cluster had instead of by
//! `pail`. This turns those into a `manifest.tbl`, after which `parse` cannot
//! tell the difference.
//!
//! Nothing is copied. The tables are found where they already are, under
//! `outputs/<pipeline>/results/`, and the only thing written is the manifest
//! beside them -- a result set large enough to be worth running elsewhere is
//! too large to want a second copy of.
//!
//! What a run was is read out of its filename, since that is already the one
//! place every pipeline in this crate records it. See [`Run::parse`] for the
//! grammar.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::Parser;

use util::manifest;
use util::tbl;
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

    /// Charge one .time file to every shard, for a search too long to have
    /// timed where it ran. Names a run, or every run without the run=
    #[arg(long, value_name = "[run=]one.time")]
    time: Vec<String>,

    /// Replace an existing manifest.tbl
    #[arg(long)]
    force: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let dir = crate::parse::pipeline(&args.pipeline)?;
    let results = dir.join("results");
    let out = dir.join("manifest.tbl");

    ensure!(
        results.is_dir(),
        "no results directory at {}; the tables go there",
        results.display()
    );

    // a pipeline that ran here wrote its own manifest, and overwriting it
    // would throw away the timings and the argv of a real run
    ensure!(
        args.force || !out.exists(),
        "{} already exists; pass --force to replace it",
        out.display()
    );

    if let Some(tool) = &args.tool {
        ensure!(
            TOOLS.contains(&tool.as_str()),
            "--tool {tool:?} is not one of {}",
            TOOLS.join(", ")
        );
    }

    let extrapolate = Extrapolate::parse(&args.time)?;

    let mut runs = collect(&results, args.tool.as_deref())?;
    ensure!(
        !runs.is_empty(),
        "no result tables in {}; expected <name>.<shard>.tbl",
        results.display()
    );

    for run in &mut runs {
        run.check_domains(&results)?;
        run.check_targets()?;
        run.read_timings(&results, extrapolate.get(&run.name))?;
    }

    write(&out, &runs)?;

    println!("wrote {}", out.display());
    for run in &runs {
        let wall: f64 = run
            .shards
            .iter()
            .filter_map(|s| s.timing.as_ref())
            .map(|t| t.wall_s)
            .sum();
        println!(
            "  {:<24} {:<7} {:>5} shard(s) {wall:>10.2}s",
            run.name,
            run.tool,
            run.shards.len()
        );
    }

    Ok(())
}

/// One run's results: what it was, and one entry per shard it covered.
struct Run {
    name: String,
    tool: String,

    /// Whatever the name said past the tool, which is what tells one run of a
    /// tool apart from another.
    params: Vec<(String, String)>,

    shards: Vec<Shard>,
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

        let tool = tool.or_else(|| fallback.map(str::to_string)).with_context(|| {
            format!("{name:?} names no tool; pass --tool to say what produced it")
        })?;

        Ok((
            Run {
                name: name.to_string(),
                tool,
                params,
                shards: Vec::new(),
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

    /// Reads each shard's `.time`, or charges every shard the one `extrapolate`
    /// names.
    ///
    /// Extrapolating is for a search whose whole point is that it was too
    /// expensive to run here: one shard is timed on a quiet machine and the
    /// rest are taken to have cost the same.
    fn read_timings(
        &mut self,
        results: &Path,
        extrapolate: Option<&Path>,
    ) -> anyhow::Result<()> {
        if let Some(path) = extrapolate {
            let shared = time::read(path)?;

            for shard in &mut self.shards {
                shard.timing = Some(Timing { ..shared });
            }

            return Ok(());
        }

        for shard in &mut self.shards {
            let path = results.join(format!("{}.{}.time", self.name, shard.name));

            ensure!(
                path.is_file(),
                "no {}\n\n\
                 every table wants a .time beside it, or --time {}=<file> to \
                 charge one to every shard",
                path.display(),
                self.name
            );

            shard.timing = Some(time::read(&path)?);
        }

        Ok(())
    }
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

/// Writes the manifest, in the shape `parse` reads back.
///
/// Assembled by hand rather than by a [`pail::Table`] sink, since no pipeline
/// ran: the commands happened somewhere else and all that came back is what
/// they produced and what they cost. Only the columns a reader uses are here.
fn write(path: &Path, runs: &[Run]) -> anyhow::Result<()> {
    // every setting any run recorded, so a pipeline that swept two knobs gets
    // two columns and one that swept none gets none
    let mut keys: Vec<&str> = runs
        .iter()
        .flat_map(|run| run.params.iter().map(|(key, _)| key.as_str()))
        .collect();
    keys.sort_unstable();
    keys.dedup();

    let mut headers = vec![
        manifest::NAME.to_string(),
        manifest::TOOL.to_string(),
        manifest::SHARD.to_string(),
    ];
    headers.extend(keys.iter().map(|key| key.to_string()));
    headers.extend(
        ["wall(s)", "user(s)", "sys(s)", "cpu(%)", "max_rss", "exit"].map(str::to_string),
    );

    let mut rows: Vec<Vec<String>> = Vec::new();

    for run in runs {
        for shard in &run.shards {
            let timing = shard
                .timing
                .as_ref()
                .expect("read_timings fills every shard");

            let mut row = vec![run.name.clone(), run.tool.clone(), shard.name.clone()];

            row.extend(keys.iter().map(|key| {
                run.params
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, value)| value.clone())
                    .unwrap_or_else(|| "-".to_string())
            }));

            row.push(format!("{:.2}", timing.wall_s));
            row.push(seconds(timing.user_s));
            row.push(seconds(timing.sys_s));
            row.push(timing.cpu_pct.map_or_else(dash, |p| format!("{p:.0}%")));
            row.push(timing.max_rss_kb.map_or_else(dash, rss));

            // a table that exists is a search that produced output, so a
            // format with nothing to say about the status is taken at that
            row.push(timing.exit.unwrap_or(0).to_string());

            rows.push(row);
        }
    }

    tbl::write(
        path,
        tbl::Table {
            meta: "",
            headers: &headers,
            rows: &rows,
            ragged_last: false,
        },
    )
}

fn dash() -> String {
    "-".to_string()
}

fn seconds(value: Option<f64>) -> String {
    value.map_or_else(dash, |s| format!("{s:.2}"))
}

/// Kilobytes as the units `pail` writes them in, so a manifest reads the same
/// however it was made.
fn rss(kb: u64) -> String {
    const MIB: f64 = 1024.0;

    match kb as f64 / MIB {
        mib if mib >= MIB => format!("{:.2}GiB", mib / MIB),
        mib => format!("{mib:.3}MiB"),
    }
}
