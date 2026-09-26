//! Turning a finished ladder into the two tables the figure and a reader want.
//!
//! `scaling.tbl` is one row per unit per arm per rung: the median over reps of
//! wall clock and core-seconds, the largest resident set, the time perfect
//! scaling from the arm's lowest rung would take, the percentage of that
//! perfect speedup the measured time reached, and the highest 1-minute load
//! seen either side of any rep.
//!
//! `agree.tbl` is one row per run: how its hits compare with the first rep of
//! its arm's lowest rung. Thread count should change nothing a tool reports,
//! so every nonzero cell there is a question for that tool.
//!
//! What ran comes out of `ledger.tbl`, so a run that failed is left out with a
//! warning rather than read as a missing file.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use clap::Parser;

use libsail::collection::Iterable;
use libsail::tbl::blast::BlastTable;
use libsail::tbl::hmmer::HmmerTable;
use libsail::tbl::nail::NailTable;
use libsail::tbl::{HitColumns, HitParser, Table};

use util::ledger::{self, Ledger, Row};
use util::manifest;

use crate::run::{ARM, Arm, LOAD, REP, THREADS};

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Where the run wrote. Defaults to what the label names
    #[arg(long, value_name = "dir")]
    out: Option<PathBuf>,
}

/// One finished run, with its place on the ladder read out of its settings.
struct Run<'a> {
    row: &'a Row,
    arm: String,
    threads: usize,
    rep: usize,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let out = args.out.unwrap_or_else(|| paths.run.clone());

    let ran = Ledger::load(&out)?;
    ledger::warn(ran.failed(), "run(s)");

    let mut runs = Vec::new();
    for row in ran.runs() {
        let param = |key: &str| {
            row.params
                .get(key)
                .with_context(|| format!("run {:?} has no {key} setting", row.name))
        };

        runs.push(Run {
            row,
            arm: param(ARM)?.clone(),
            threads: param(THREADS)?.parse().context("threads")?,
            rep: param(REP)?.parse().context("rep")?,
        });
    }
    if runs.is_empty() {
        bail!("no finished runs in {}", out.display());
    }

    // (shard, arm) then threads, so each arm's ladder comes out ascending
    let mut ladders: BTreeMap<(&str, &str), BTreeMap<usize, Vec<&Run>>> = BTreeMap::new();
    for run in &runs {
        ladders
            .entry((&run.row.shard, &run.arm))
            .or_default()
            .entry(run.threads)
            .or_default()
            .push(run);
    }

    // hmmer-split at one thread is one part at --cpu 1, which is the
    // command hmmer's own one-thread rung ran, so it takes that point as
    // its first rather than starting its ladder at 2
    let (hmmer, split) = (Arm::Hmmer.name(), Arm::HmmerSplit.name());
    let ones: Vec<(&str, Vec<&Run>)> = ladders
        .iter()
        .filter(|((_, arm), _)| *arm == hmmer)
        .filter_map(|(&(shard, _), rungs)| rungs.get(&1).map(|one| (shard, one.clone())))
        .collect();
    for (shard, one) in ones {
        if let Some(rungs) = ladders.get_mut(&(shard, split)) {
            rungs.entry(1).or_insert(one);
        }
    }

    std::fs::create_dir_all(&paths.analysis)
        .with_context(|| format!("failed to create {}", paths.analysis.display()))?;

    let load = load(&out.join(LOAD))?;
    let scaling = paths.analysis.join("scaling.tbl");
    write_scaling(&scaling, &ladders, &load)?;
    println!("wrote {}", scaling.display());

    let agree = paths.analysis.join("agree.tbl");
    write_agree(&agree, &out.join("results"), &ladders)?;
    println!("wrote {}", agree.display());

    Ok(())
}

fn median(mut xs: Vec<f64>) -> Option<f64> {
    xs.sort_by(f64::total_cmp);
    let n = xs.len();
    match n {
        0 => None,
        _ if n % 2 == 1 => Some(xs[n / 2]),
        _ => Some((xs[n / 2 - 1] + xs[n / 2]) / 2.0),
    }
}

/// The highest 1-minute load read against each (run, shard).
///
/// A run from before the load log existed has none, and reads as `-`.
fn load(path: &Path) -> anyhow::Result<HashMap<(String, String), f64>> {
    let mut out: HashMap<(String, String), f64> = HashMap::new();
    if !path.is_file() {
        return Ok(out);
    }

    let table =
        toil::Table::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    let at = |label: &str| {
        table
            .index(label)
            .with_context(|| format!("{} has no {label} column", path.display()))
    };
    let (name, shard, load1) = (at("name")?, at("shard")?, at("load1")?);

    for cells in table.rows() {
        let cell = |i: usize| cells.get(i).unwrap_or_default();
        let value: f64 = cell(load1)
            .parse()
            .with_context(|| format!("{:?} in {} is not a load", cell(load1), path.display()))?;

        let max = out
            .entry((cell(name).to_string(), cell(shard).to_string()))
            .or_insert(value);
        *max = max.max(value);
    }

    Ok(out)
}

type Ladders<'a> = BTreeMap<(&'a str, &'a str), BTreeMap<usize, Vec<&'a Run<'a>>>>;

fn write_scaling(
    path: &Path,
    ladders: &Ladders,
    load: &HashMap<(String, String), f64>,
) -> anyhow::Result<()> {
    let mut table = toil::Table::new(toil::Schema::new([
        "unit",
        "arm",
        "threads",
        "reps",
        "wall_s",
        "cpu_s",
        "max_rss",
        "ideal_s",
        "pct_ideal",
        "load",
    ]));
    table.meta("wall_s and cpu_s are medians over reps");
    table.meta(
        "ideal_s is the arm's lowest rung scaled perfectly: wall_s x lowest threads / threads",
    );

    let dash = || "-".to_string();

    for ((shard, arm), rungs) in ladders {
        let wall = |runs: &[&Run]| median(runs.iter().filter_map(|r| r.row.wall_s).collect());
        let (base_threads, base_wall) = rungs
            .iter()
            .next()
            .map(|(&n, runs)| (n, wall(runs)))
            .expect("a ladder holds at least the rung that made it");

        for (&threads, runs) in rungs {
            let wall_s = wall(runs);
            let cpu_s = median(runs.iter().filter_map(|r| r.row.cpu_s).collect());
            let rss = runs.iter().filter_map(|r| r.row.max_rss_kb).max();
            let ideal = base_wall.map(|b| b * base_threads as f64 / threads as f64);
            let pct = ideal.zip(wall_s).map(|(i, w)| 100.0 * i / w);
            let max_load = runs
                .iter()
                .filter_map(|r| load.get(&(r.row.name.clone(), r.row.shard.clone())))
                .copied()
                .reduce(f64::max);

            table.row([
                shard.to_string(),
                arm.to_string(),
                threads.to_string(),
                runs.len().to_string(),
                wall_s.map_or_else(dash, |x| format!("{x:.2}")),
                cpu_s.map_or_else(dash, |x| format!("{x:.2}")),
                rss.map_or_else(dash, |x| x.to_string()),
                ideal.map_or_else(dash, |x| format!("{x:.2}")),
                pct.map_or_else(dash, |x| format!("{x:.1}")),
                max_load.map_or_else(dash, |x| format!("{x:.2}")),
            ]);
        }
    }

    table
        .write(path)
        .with_context(|| format!("failed to write {}", path.display()))
}

/// Each (query, target) a table reports, and the best score it gave the pair.
fn hits(tool: &str, path: &Path) -> anyhow::Result<HashMap<(String, String), f32>> {
    fn read<C: HitColumns>(path: &Path) -> anyhow::Result<HashMap<(String, String), f32>> {
        let table = Table::<HitParser<C>>::open(path)
            .with_context(|| format!("failed to open {}", path.display()))?;

        let mut out: HashMap<(String, String), f32> = HashMap::new();
        for hit in table.iter() {
            // mmseqs reports a pair once per alignment, so a pair can repeat
            let best = out
                .entry((hit.query.clone(), hit.target.clone()))
                .or_insert(hit.score);
            *best = best.max(hit.score);
        }
        Ok(out)
    }

    match tool {
        "nail" => read::<NailTable>(path),
        "mmseqs" => read::<BlastTable>(path),
        "hmmer" => read::<HmmerTable>(path),
        other => bail!("no reader for {other}'s tables"),
    }
}

fn write_agree(path: &Path, results: &Path, ladders: &Ladders) -> anyhow::Result<()> {
    let mut table = toil::Table::new(toil::Schema::new([
        "unit",
        "arm",
        "threads",
        "rep",
        "hits",
        "missing",
        "extra",
        "rescored",
        "max_delta",
    ]));
    table.meta("against rep 1 of each arm's lowest rung, or its lowest rep there");

    for ((shard, arm), rungs) in ladders {
        let base = rungs
            .values()
            .next()
            .and_then(|runs| runs.iter().min_by_key(|r| r.rep))
            .expect("a ladder holds at least the run that made it");
        let table_of = |run: &Run| manifest::table_path(results, &run.row.name, &run.row.shard);
        let base_hits = hits(&base.row.tool, &table_of(base))?;

        for runs in rungs.values() {
            let mut runs = runs.clone();
            runs.sort_by_key(|r| r.rep);

            for run in runs {
                let found = hits(&run.row.tool, &table_of(run))?;

                let missing = base_hits.keys().filter(|k| !found.contains_key(*k)).count();
                let extra = found.keys().filter(|k| !base_hits.contains_key(*k)).count();
                let deltas: Vec<f32> = found
                    .iter()
                    .filter_map(|(k, s)| base_hits.get(k).map(|b| (s - b).abs()))
                    .collect();
                let rescored = deltas.iter().filter(|&&d| d > 0.0).count();
                let max_delta = deltas.iter().copied().fold(0.0f32, f32::max);

                table.row([
                    shard.to_string(),
                    arm.to_string(),
                    run.threads.to_string(),
                    run.rep.to_string(),
                    found.len().to_string(),
                    missing.to_string(),
                    extra.to_string(),
                    rescored.to_string(),
                    format!("{max_delta:.2}"),
                ]);
            }
        }
    }

    table
        .write(path)
        .with_context(|| format!("failed to write {}", path.display()))
}
