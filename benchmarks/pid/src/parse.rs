//! Turning a finished run into the tables the plot scripts consume.
//!
//! What ran comes out of `manifest.tbl` -- the run's name, which tool produced
//! its table, which query it searched, and what it cost -- rather than out of
//! the results directory's filenames. That is what lets a run be renamed, or a
//! tool added, without this file learning about it.
//!
//! Truth here is in the benchmark itself: `truth.tbl` says which pair is
//! which and at what identity, and a target named `decoy…` is one. There is no
//! calibration and no reference tool. A decoy whose original the query's
//! family hits is the reversal of a homolog, and recall skips it; `reject` is
//! what finds those, with hmmsearch judging for every tool.

use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
};

use libsail::collection::Iterable;
use libsail::seq::fasta::IndexedFasta;
use libsail::seq::p7hmm::Hmm;
use libsail::tbl::blast::BlastTable;
use libsail::tbl::hmmer::HmmerTable;
use libsail::tbl::nail::NailTable;
use libsail::tbl::{Hit, HitColumns, HitParser, Table};

use anyhow::Context;
use clap::{Parser, Subcommand};
use util::ledger::{self, Ledger};
use util::manifest;

use crate::run::MODE;

const PRECISION: usize = 4;
const FIXED_FPR: f32 = 0.01;

/// The E-value at or under which a family's hit against a decoy's original
/// makes that decoy a reject for that family, as cutoffs has it.
const REJECT_E: f64 = 1e-3;

/// The E-value a run has to beat to sit at [`FIXED_FPR`], which is the score
/// of the decoy that many places down its own list.
//
// TODO: why decoy_cnt + 1 rather than decoy_cnt?
fn threshold(tbl: &HitTable2, positive_cnt: usize) -> f64 {
    let decoy_cnt = (positive_cnt as f32 * FIXED_FPR).ceil() as usize;

    match tbl.adjusted_decoys.get(decoy_cnt + 1) {
        Some(e_value) => *e_value,
        None => {
            println!(
                "warning: not enough decoys to produce E-value threshold for: {}",
                tbl.name
            );
            f64::INFINITY
        }
    }
}

fn e_value_cmp(a: f64, b: f64) -> std::cmp::Ordering {
    a.partial_cmp(&b).expect("NaN encountered in E-value cmp")
}

/// Where an analysis writes. Every one of them takes it.
#[derive(Parser)]
pub struct Which {
    /// Where the tables go. Defaults to what the label names
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,
}

impl Which {
    fn out_dir(&self, analysis: &Path) -> PathBuf {
        self.out.clone().unwrap_or_else(|| analysis.to_path_buf())
    }
}

#[derive(Parser)]
pub struct RecallArgs {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    #[command(flatten)]
    which: Which,
}

#[derive(Parser)]
pub struct CellsArgs {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    #[command(flatten)]
    which: Which,

    /// Which run's table to read cell fractions from. Only nail reports them
    #[arg(long, value_name = "NAME", default_value = "nail-s12.0-ms2000.prf")]
    run: String,
}

#[derive(Parser)]
pub struct ScoreArgs {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Two nail tables to correlate, by run name. There is no --full-dp run in
    /// the sweep today, so these are given rather than assumed
    #[arg(long, value_name = "NAME")]
    full: String,

    #[arg(long, value_name = "NAME")]
    sparse: String,

    #[command(flatten)]
    which: Which,
}

#[derive(Subcommand)]
pub enum Cmd {
    Recall(RecallArgs),
    Cells(CellsArgs),
    Score(ScoreArgs),
}

pub fn main(cmd: Cmd, paths: &crate::Paths) -> anyhow::Result<()> {
    match cmd {
        Cmd::Recall(args) => {
            recall(args, paths)?;
        }
        Cmd::Cells(args) => {
            cells(args, paths)?;
        }
        Cmd::Score(args) => {
            score(args, paths)?;
        }
    }

    Ok(())
}

fn score(args: ScoreArgs, paths: &crate::Paths) -> anyhow::Result<()> {
    let results = paths.search().join("results");

    let read = |name: &str| -> anyhow::Result<HashMap<(String, String), f32>> {
        let path = manifest::table_path(&results, name, "");
        let tbl = Table::<HitParser<NailTable>>::open(&path)
            .with_context(|| format!("failed to open {}", path.display()))?;

        // one entry per pair; a repeated pair keeps its last row
        Ok(tbl
            .iter()
            .map(|h| ((h.query.clone(), h.target.clone()), h.score))
            .collect())
    };

    let full_tbl = read(&args.full)?;
    let sparse_tbl = read(&args.sparse)?;

    let intersection = full_tbl
        .keys()
        .filter(|k| sparse_tbl.contains_key(*k))
        .collect::<Vec<_>>();

    let figures = args.which.out_dir(&paths.analysis);
    std::fs::create_dir_all(&figures)?;

    let mut table = toil::Table::new(toil::Schema::new([
        toil::Column::new("full").fixed(1),
        toil::Column::new("sparse").fixed(1),
    ]));
    for k in intersection {
        let x = full_tbl.get(k).expect("present by intersection");
        let y = sparse_tbl.get(k).expect("present by intersection");
        table.row([*x, *y]);
    }

    let path = figures.join("score.tbl");
    table
        .write(&path)
        .with_context(|| format!("failed to write {}", path.display()))
}

fn cells(args: CellsArgs, paths: &crate::Paths) -> anyhow::Result<()> {
    let table = manifest::table_path(&paths.search().join("results"), &args.run, "");

    let hits = libsail::tbl::NailRows::open(&table)
        .with_context(|| format!("failed to read {}", table.display()))?;

    // read out of the files rather than shelled out to hmmstat and
    // esl-seqstat: neither is a dependency this benchmark declares, and both
    // were being found on PATH rather than through `tools`
    let inp = crate::Inputs::open(&paths.set)?;
    let query_lens: HashMap<String, usize> = Hmm::open(&inp.query_hmm)
        .with_context(|| format!("failed to parse {}", inp.query_hmm.display()))?
        .iter()
        .map(|model| {
            (
                String::from_utf8_lossy(&model.header.name).into_owned(),
                model.header.leng,
            )
        })
        .collect();

    let mut target_lens: HashMap<String, usize> = HashMap::new();
    let target_fa = IndexedFasta::open(&inp.target_fa)
        .with_context(|| format!("failed to open {}", inp.target_fa.display()))?;
    for rec in target_fa.iter() {
        target_lens.insert(rec.name_str()?.to_string(), rec.seq.len());
    }

    let figures = args.which.out_dir(&paths.analysis);
    std::fs::create_dir_all(&figures)?;

    let points = || toil::Table::new(toil::Schema::new(["cells", "fraction"]));
    let (mut true_out, mut decoy_out) = (points(), points());

    hits.iter().try_for_each(|h| -> anyhow::Result<()> {
        let intended_query = h
            .target
            .split('|')
            .next()
            .with_context(|| format!("failed to split query from: {}", h.target))?;

        let qlen = query_lens
            .get(&h.query)
            .with_context(|| format!("no query len for: {}", h.query))?;

        let tlen = target_lens
            .get(&h.target)
            .with_context(|| format!("no target len for: {}", h.target))?;

        let x = toil::Cell::from(qlen * tlen);
        let y = toil::Cell::from(h.cell_fraction);

        if h.target.starts_with("decoy") {
            decoy_out.row([x, y]);
        } else if h.query == intended_query {
            true_out.row([x, y]);
        }

        Ok(())
    })?;

    for (table, name) in [(true_out, "cells.true.tbl"), (decoy_out, "cells.decoy.tbl")] {
        let path = figures.join(name);
        table
            .write(&path)
            .with_context(|| format!("failed to write {}", path.display()))?;
    }

    Ok(())
}

fn recall(args: RecallArgs, paths: &crate::Paths) -> anyhow::Result<()> {
    let start = std::time::Instant::now();

    let inp = crate::Inputs::open(&paths.set)?;
    let benchmark = Benchmark::new(inp.truth)?;
    let rejected = rejections(&paths.reject())?;
    let data = RecallData::new(&paths.search(), &benchmark, &rejected)?;

    let figures = args.which.out_dir(&paths.analysis);
    std::fs::create_dir_all(&figures)?;

    data.write_roc(&figures.join("roc.tbl"))?;
    data.write_pid(&figures.join("pid.tbl"))?;
    data.write_runtime(&figures.join("time.tbl"))?;

    println!("recall data took: {:?}", start.elapsed());
    Ok(())
}

pub struct BenchmarkEntry {
    pid: usize,
    query: String,
    family: String,
}

pub struct Benchmark {
    entries: Vec<BenchmarkEntry>,
    idx_by_target: HashMap<String, usize>,
}

impl Benchmark {
    pub fn new<P: AsRef<Path>>(tbl_path: P) -> anyhow::Result<Self> {
        let path = tbl_path.as_ref();
        let table = toil::Table::read(path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        let at = |label: &str| {
            table
                .index(label)
                .with_context(|| format!("{} has no {label} column", path.display()))
        };
        let (pid_at, family_at, target_at, query_at) =
            (at("identity")?, at("family")?, at("target")?, at("query")?);

        let mut entries = vec![];
        let mut idx_by_target: HashMap<String, usize> = HashMap::new();

        for cells in table.rows() {
            let cell = |i: usize| {
                cells
                    .get(i)
                    .with_context(|| format!("{} has a row with an empty cell", path.display()))
            };

            let pid = cell(pid_at)?
                .parse::<usize>()
                .context("identity is not a whole percentage")?;
            let family = cell(family_at)?.to_string();
            let target = cell(target_at)?.to_string();
            let query = cell(query_at)?.to_string();

            let entry = BenchmarkEntry {
                pid,
                query: query.clone(),
                family: family.clone(),
            };

            idx_by_target.insert(target, entries.len());
            entries.push(entry);
        }

        Ok(Self {
            entries,
            idx_by_target,
        })
    }
}

impl Benchmark {
    fn entry_by_target(&self, target: &str) -> &BenchmarkEntry {
        &self.entries[*self.idx_by_target.get(target).unwrap()]
    }
}

/// One finished run: what it was called, how to read its table, which query it
/// searched, and what it cost.
pub struct Run {
    pub name: String,
    pub tool: String,
    /// `prf`, `cons` or `seq`, as the run recorded it.
    pub mode: String,
    wall_s: f32,
    table: PathBuf,
}

impl Run {
    /// The hits it reported. Which reader to use comes off the `tool` field --
    /// mmseqs, last, blast and diamond all write blast's tabular format,
    /// whatever wrote it.
    pub fn hits(&self) -> anyhow::Result<Vec<Hit>> {
        match self.tool.as_str() {
            "nail" => read_hits::<NailTable>(&self.table),
            "hmmer" | "phmmer" => read_hits::<HmmerTable>(&self.table),
            _ => read_hits::<BlastTable>(&self.table),
        }
        .with_context(|| format!("failed to read {}", self.table.display()))
    }
}

/// Every row of a table read as layout `C`.
fn read_hits<C: HitColumns>(path: &Path) -> libsail::Result<Vec<Hit>> {
    Ok(Table::<HitParser<C>>::open(path)?.iter().cloned().collect())
}

/// The runs a pipeline finished, in the order it declared them.
///
/// Read out of `ledger.tbl` rather than by globbing the results directory,
/// which
/// is what keeps the runs table itself from looking like a hit table and
/// what makes a run's tool and mode facts rather than guesses.
///
/// Several commands can share a name -- mmseqs' search and its conversion are
/// one run, and psiblast's per-family calls are one run done a family at a
/// time -- and they arrive here already folded into one row.
pub fn runs(dir: &Path) -> anyhow::Result<Vec<Run>> {
    let ran = Ledger::load(dir)?;
    ledger::warn(ran.failed(), "command(s)");

    let results = dir.join("results");
    let out: Vec<Run> = ran
        .columns()?
        .into_iter()
        .map(|column| {
            let mode = column
                .params
                .get(MODE)
                .with_context(|| format!("run {:?} has no mode", column.name))?;

            anyhow::ensure!(
                matches!(mode.as_str(), "prf" | "cons" | "seq"),
                "unknown search mode {mode:?} in ledger.tbl"
            );

            Ok(Run {
                table: manifest::table_path(&results, &column.name, ""),
                mode: mode.clone(),
                name: column.name,
                tool: column.tool,
                wall_s: column.wall_s as f32,
            })
        })
        .collect::<anyhow::Result<_>>()?;

    Ok(out)
}

/// Every (family, decoy) pair whose original the judge found, which is every
/// pair a run's decoy hit is dropped for.
///
/// A pair is settled for the one family whose profile hit the original. A
/// sequence query is judged by its family's profile. Nothing here pools a
/// decoy across families.
pub fn rejections(dir: &Path) -> anyhow::Result<Rejections> {
    let tested = toil::Table::read(dir.join(TESTED))
        .with_context(|| format!("failed to read {}", dir.join(TESTED).display()))
        .with_context(|| format!("no rejections in {}; run `pid reject` first", dir.display()))?;
    let column = |label: &str| {
        tested
            .index(label)
            .with_context(|| format!("{} has no {label} column", dir.join(TESTED).display()))
    };
    let (tool, decoys) = (column("tool")?, column("decoys")?);

    let settled = tested
        .rows()
        .iter()
        .find(|row| row.get(tool) == Some("all"))
        .with_context(|| format!("{} has no total row", dir.join(TESTED).display()))?;

    if settled.get(decoys) == Some("0") {
        return Ok(HashSet::new());
    }

    // decoys to settle and no finished judge would read as nothing to reject,
    // which is the answer that inflates every run's false positives
    let judge = runs(dir)?
        .into_iter()
        .find(|run| run.name == crate::reject::JUDGE)
        .with_context(|| format!("no finished judge in {}", dir.display()))?;

    Ok(judge
        .hits()?
        .into_iter()
        .filter(|h| h.e_value <= REJECT_E)
        .map(|h| (h.query, h.target))
        .collect())
}

/// The rejected (family, decoy) pairs.
pub type Rejections = HashSet<(String, String)>;

/// What `reject` writes before it searches: how many decoys each tool had to
/// settle, including the ones with none, and a total row for the union.
pub const TESTED: &str = "tested.tbl";

/// A run's hits reduced to the best E-value per (target, query).
fn best(hits: &[Hit]) -> HashMap<(&str, &str), f64> {
    // a run can report one pair more than once; the best of them is the
    // one a threshold would see
    let mut best: HashMap<(&str, &str), f64> = HashMap::new();
    for hit in hits {
        best.entry((&hit.target, &hit.query))
            .and_modify(|e| *e = e.min(hit.e_value))
            .or_insert(hit.e_value);
    }
    best
}

/// A query's family, and the id a decoy list is keyed by.
fn query_parts(query: &str) -> (&str, &str) {
    // a query is <family> for a profile and <family>|<id> for a
    // sequence, and the second is what a decoy list is keyed by
    match query.split_once('|') {
        Some((family, id)) => (family, id),
        None => (query, query),
    }
}

/// The identity a true pair was drawn at, if `target` is a true target of
/// `query`.
fn truth_of(bm: &Benchmark, target: &str, query: &str, search_type: &str) -> Option<usize> {
    let (family, id) = query_parts(query);

    // a target's name carries the pair it is: <family>|<id>|<pid>%
    let target = target
        .split('|')
        .nth(1)
        .unwrap_or_else(|| panic!("target {target:?} names no sequence"));

    let entry = bm.entry_by_target(target);

    let matched = match search_type {
        "seq" => id == entry.query,
        _ => family == entry.family,
    };

    matched.then_some(entry.pid)
}

/// The decoys a run ranked at or above its worst true pair.
///
/// These are the only decoys that can move a point of its ROC: below the
/// last true pair the curve has reached its full recall, and removing a decoy
/// there only shortens the flat tail.
pub fn tested<'a>(hits: &'a [Hit], bm: &Benchmark, search_type: &str) -> HashSet<&'a str> {
    let best = best(hits);

    let floor = best
        .iter()
        .filter(|((target, query), _)| {
            !target.starts_with("decoy") && truth_of(bm, target, query, search_type).is_some()
        })
        .map(|(_, e_value)| *e_value)
        .max_by(|a, b| e_value_cmp(*a, *b));

    let Some(floor) = floor else {
        return HashSet::new();
    };

    // at the floor rather than under it: the ROC counts a true pair only once
    // its E-value is strictly under a decoy's, so a tie outranks it
    best.into_iter()
        .filter(|((target, _), e_value)| target.starts_with("decoy") && *e_value <= floor)
        .map(|((target, _), _)| target)
        .collect()
}

struct HitTable2 {
    name: String,
    /// Every true pair this run reported, as the identity it was drawn at and
    /// the E-value the run gave it, worst last.
    positives: Vec<(usize, f64)>,
    /// The decoys this run's own queries reported, in the benchmark's order.
    ///
    /// What the ROC walks: a false positive only counts against a query that
    /// was actually asked, so the list is rebuilt per benchmark entry rather
    /// than taken as everything the run called a decoy.
    adjusted_decoys: Vec<f64>,
}

impl HitTable2 {
    fn new(
        hits: &[Hit],
        name: &str,
        bm: &Benchmark,
        search_type: &str,
        rejected: &Rejections,
    ) -> Self {
        let mut positives: Vec<(usize, f64)> = vec![];
        let mut decoys_by_query: HashMap<&str, Vec<f64>> = HashMap::new();

        for ((target, query), e_value) in best(hits) {
            if target.starts_with("decoy") {
                let family = query_parts(query).0;
                if rejected.contains(&(family.to_string(), target.to_string())) {
                    continue;
                }

                decoys_by_query
                    .entry(query_parts(query).1)
                    .or_default()
                    .push(e_value);
                continue;
            }

            if let Some(pid) = truth_of(bm, target, query, search_type) {
                positives.push((pid, e_value));
            }
        }

        // a false positive only counts against a query that was asked, so the
        // decoys are gathered per benchmark entry rather than taken whole
        let mut adjusted_decoys: Vec<f64> = bm
            .entries
            .iter()
            .flat_map(|entry| {
                let asked = match search_type {
                    "cons" => &format!("{}-consensus", entry.family),
                    "seq" => &entry.query,
                    _ => &entry.family,
                };

                decoys_by_query
                    .get(asked.as_str())
                    .cloned()
                    .unwrap_or_default()
            })
            .collect();

        positives.sort_by(|a, b| e_value_cmp(a.1, b.1));
        adjusted_decoys.sort_by(|a, b| e_value_cmp(*a, *b));

        Self {
            name: name.to_string(),
            positives,
            adjusted_decoys,
        }
    }
}

struct RecallData {
    tables: Vec<HitTable2>,
    times: Vec<f32>,
    bin_sizes: Vec<usize>,
    positive_cnt: usize,
}

impl RecallData {
    fn new(dir: &Path, bm: &Benchmark, rejected: &Rejections) -> anyhow::Result<Self> {
        let mut pid_bin_tot_cnts = vec![];
        bm.entries.iter().map(|e| e.pid).for_each(|pid| {
            if pid >= pid_bin_tot_cnts.len() {
                pid_bin_tot_cnts.resize(pid + 1, 0);
            }
            pid_bin_tot_cnts[pid] += 1;
        });

        let mut tables = vec![];
        let mut times = vec![];

        let ran = runs(dir)?;
        anyhow::ensure!(!ran.is_empty(), "no finished runs in {}", dir.display());

        for run in ran {
            tables.push(HitTable2::new(
                &run.hits()?,
                &run.name,
                bm,
                &run.mode,
                rejected,
            ));
            times.push(run.wall_s);
        }

        Ok(Self {
            tables,
            times,
            bin_sizes: pid_bin_tot_cnts,
            positive_cnt: bm.entries.len(),
        })
    }

    /// One row per run per identity bin: what fraction of that bin's true
    /// pairs the run kept at [`FIXED_FPR`].
    fn write_pid(&self, path: &Path) -> anyhow::Result<()> {
        let mut rows = vec![];

        for tbl in &self.tables {
            let cutoff = threshold(tbl, self.positive_cnt);

            let mut kept = vec![0usize; self.bin_sizes.len()];
            for (pid, e_value) in &tbl.positives {
                if *e_value <= cutoff {
                    kept[*pid] += 1;
                }
            }

            for (pid, (&kept, &size)) in kept.iter().zip(&self.bin_sizes).enumerate() {
                if size == 0 {
                    continue;
                }

                assert!(
                    kept <= size,
                    "bin count > bin size: {pid}% | {kept} > {size}"
                );

                rows.push(vec![
                    tbl.name.clone(),
                    pid.to_string(),
                    size.to_string(),
                    format!("{:.p$}", kept as f64 / size as f64, p = PRECISION),
                ]);
            }
        }

        let mut table = toil::Table::new(toil::Schema::new(["run", "pid", "n", "recall"]));
        table.meta("fpr", [format!("{FIXED_FPR}")]);
        for row in rows {
            table.row(row);
        }

        table
            .write(path)
            .with_context(|| format!("failed to write {}", path.display()))
    }

    /// One row per run per change point of its ROC curve.
    fn write_roc(&self, path: &Path) -> anyhow::Result<()> {
        let p = 10.0f64.powi(PRECISION as i32);
        let mut rows = vec![];

        for tbl in &self.tables {
            let mut e_values = tbl.adjusted_decoys.clone();
            e_values.push(f64::INFINITY);

            // each decoy admitted is one more false positive, and what moves
            // is how many true pairs came in under it
            let mut counts = vec![];
            let mut found = 0usize;
            for e in e_values {
                found += tbl.positives[found..]
                    .iter()
                    .take_while(|(_, e_value)| *e_value < e)
                    .count();

                counts.push(found);
            }

            let all: Vec<(f64, f64)> = counts
                .into_iter()
                .enumerate()
                .map(|(i, c)| {
                    let x = i as f64 / self.positive_cnt as f64;
                    let y = c as f64 / self.positive_cnt as f64;
                    ((x * p).round() / p, (y * p).round() / p)
                })
                .collect();

            // only where the curve turns: a run of points at one recall draws
            // the same line as its two ends
            let mut points = vec![all[0]];
            for pair in all.windows(2) {
                if pair[0].1 != pair[1].1 {
                    points.push(pair[0]);
                    points.push(pair[1]);
                }
            }
            points.dedup();

            for (x, y) in points {
                rows.push(vec![
                    tbl.name.clone(),
                    format!("{x:.p$}", p = PRECISION),
                    format!("{y:.p$}", p = PRECISION),
                ]);
            }
        }

        let mut table = toil::Table::new(toil::Schema::new(["run", "fpr", "recall"]));
        for row in rows {
            table.row(row);
        }

        table
            .write(path)
            .with_context(|| format!("failed to write {}", path.display()))
    }

    /// One row per run: what it cost, and what it found for the cost.
    fn write_runtime(&self, path: &Path) -> anyhow::Result<()> {
        let rows: Vec<Vec<String>> = self
            .tables
            .iter()
            .zip(&self.times)
            .map(|(tbl, wall_s)| {
                let cutoff = threshold(tbl, self.positive_cnt);

                let found = tbl
                    .positives
                    .iter()
                    .take_while(|(_, e_value)| *e_value < cutoff)
                    .count();

                vec![
                    tbl.name.clone(),
                    format!("{wall_s:.4}"),
                    format!("{:.4}", found as f64 / self.positive_cnt as f64),
                ]
            })
            .collect();

        let mut table = toil::Table::new(toil::Schema::new(["run", "wall_s", "recall"]));
        table.meta("fpr", [format!("{FIXED_FPR}")]);
        for row in rows {
            table.row(row);
        }

        table
            .write(path)
            .with_context(|| format!("failed to write {}", path.display()))
    }
}
