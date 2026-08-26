//! Turning a finished run into the tables the plot scripts consume.
//!
//! What ran comes out of `manifest.tbl` -- the run's name, which tool produced
//! its table, which query it searched, and what it cost -- rather than out of
//! the results directory's filenames. That is what lets a run be renamed, or a
//! tool added, without this file learning about it.
//!
//! Truth here is in the benchmark itself: `benchmark.tbl` says which pair is
//! which and at what identity, and a target named `decoy…` is one. There is no
//! calibration and no reference tool.

use std::{
    collections::HashMap,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use bioio::tbl::{BlastTable, Hit, HitTable, HmmerTable, NailTable};

use anyhow::{Context, bail};
use bench::manifest::{self, Manifest, Wall};
use clap::{Parser, Subcommand};

use crate::inputs::Inputs;
use crate::search::MODE;

const PRECISION: usize = 4;
const FIXED_FPR: f32 = 0.01;

fn e_value_cmp(a: &Hit2, b: &Hit2) -> std::cmp::Ordering {
    a.e_value
        .partial_cmp(&b.e_value)
        .expect("NaN encountered in E-value cmp")
}

/// Which input set was searched, and so where everything is. Every analysis
/// takes it, the way `build` and `run` do.
#[derive(Parser)]
pub struct Which {
    #[arg(short, long, default_value = "toy")]
    size: String,

    /// Where the tables go. Defaults to figures/ beside the run.
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,
}

impl Which {
    fn set(&self) -> Inputs {
        Inputs::new(&self.size)
    }

    fn out_dir(&self) -> PathBuf {
        self.out
            .clone()
            .unwrap_or_else(|| self.set().run_dir().join("figures"))
    }
}

#[derive(Parser)]
pub struct RecallArgs {
    #[command(flatten)]
    which: Which,
}

#[derive(Parser)]
pub struct CellsArgs {
    #[command(flatten)]
    which: Which,

    /// Which run's table to read cell fractions from. Only nail reports them.
    #[arg(long, value_name = "NAME", default_value = "nail-s12.0-ms2000.prf")]
    run: String,
}

#[derive(Parser)]
pub struct ScoreArgs {
    /// Two nail tables to correlate, by run name. There is no --full-dp run in
    /// the sweep today, so these are given rather than assumed.
    #[arg(long, value_name = "NAME")]
    full: String,

    #[arg(long, value_name = "NAME")]
    sparse: String,

    #[command(flatten)]
    which: Which,
}

#[derive(Parser)]
pub struct TableArgs {
    #[command(flatten)]
    which: Which,

    #[arg(short = 'e', long, default_value_t = false)]
    e_value: bool,

    #[arg(long, default_value_t = 6, default_value_if("e_value", "true", "9"))]
    min_width: usize,
}

/// Analysis subcommands for this benchmark. Kept here rather than in the `run`
/// library because what counts as a result is benchmark-specific.
#[derive(Subcommand)]
pub enum Cmd {
    Recall(RecallArgs),
    Cells(CellsArgs),
    Score(ScoreArgs),
    Table(TableArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Recall(args) => {
            recall(args)?;
        }
        Cmd::Cells(args) => {
            cells(args)?;
        }
        Cmd::Score(args) => {
            score(args)?;
        }
        Cmd::Table(args) => {
            table(args)?;
        }
    }

    Ok(())
}

fn table(args: TableArgs) -> anyhow::Result<()> {
    let set = args.which.set();
    let bm = Benchmark::new(set.benchmark_tbl())?;

    let out_dir = args.which.out_dir();
    std::fs::create_dir_all(&out_dir)?;
    let mut out = BufWriter::new(File::create(out_dir.join("results.tbl"))?);

    let mut tuples = vec![];

    for run in runs(&set.run_dir())? {
        fn true_hit_filter(hit: &Hit) -> bool {
            if hit.target.starts_with("decoy") {
                return false;
            }

            let q = hit.query.split('|').next().unwrap_or_default();
            let t = hit.target.split('|').next().unwrap_or_default();

            q == t
        }

        let hits: Vec<Hit> = run
            .hits()?
            .hits
            .into_iter()
            .filter(true_hit_filter)
            .collect();

        // a target's name carries which pair it is; everything past that is
        // the family it came from and the identity it was drawn at
        let target_of = |hit: &Hit| -> anyhow::Result<String> {
            hit.target
                .split('|')
                .nth(1)
                .map(str::to_string)
                .with_context(|| format!("target {:?} names no sequence", hit.target))
        };

        let value = |hit: &Hit| match args.e_value {
            true => hit.e_value,
            false => hit.score as f64,
        };

        let mut map: HashMap<String, f64> = HashMap::new();
        for hit in &hits {
            let target = target_of(hit)?;
            let v = value(hit);

            match run.mode {
                // one profile per family, so a pair is reported once
                SearchType::Profile | SearchType::Consensus => {
                    map.insert(target, v);
                }
                // several query sequences can reach the same target, so the
                // best of them is the one a threshold would see
                SearchType::Sequence => {
                    map.entry(target)
                        .and_modify(|best| {
                            *best = match args.e_value {
                                true => best.min(v),
                                false => best.max(v),
                            }
                        })
                        .or_insert(v);
                }
            }
        }

        // the run name is <prefix>.<mode>, and the prefix's `-` pieces are what
        // the header stacks up
        let prefix = run
            .name
            .rsplit_once('.')
            .map_or(run.name.as_str(), |(p, _)| p);
        let prefix_tokens: Vec<String> = prefix.split('-').map(str::to_string).collect();

        tuples.push((prefix_tokens, run.mode.to_string(), map));
    }

    fn cmp_component(a: &str, b: &str) -> std::cmp::Ordering {
        let na = a.chars().find(|c| c.is_ascii_digit()).map(|_| {
            a.chars()
                .skip_while(|c| !c.is_ascii_digit())
                .collect::<String>()
                .parse::<f64>()
                .unwrap()
        });

        let nb = b.chars().find(|c| c.is_ascii_digit()).map(|_| {
            b.chars()
                .skip_while(|c| !c.is_ascii_digit())
                .collect::<String>()
                .parse::<f64>()
                .unwrap()
        });

        match (na, nb) {
            (Some(a), Some(b)) => a.partial_cmp(&b).unwrap(),
            _ => a.cmp(b),
        }
    }

    fn cmp_keys(a: &[String], b: &[String]) -> std::cmp::Ordering {
        for (x, y) in a.iter().zip(b) {
            let ord = cmp_component(x, y);
            if ord != std::cmp::Ordering::Equal {
                return ord;
            }
        }

        a.len().cmp(&b.len())
    }

    tuples.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| cmp_keys(&a.0, &b.0)));

    let n_rows = tuples.iter().map(|(p, _, _)| p.len()).max().unwrap();

    let target_width = bm.entries.iter().map(|e| e.target.len()).max().unwrap();
    let fam_width = bm.entries.iter().map(|e| e.family.len()).max().unwrap();

    let widths = tuples
        .iter()
        .map(|(p, _, _)| p.iter().map(|s| s.len().max(args.min_width)).max().unwrap())
        .collect::<Vec<_>>();

    let total_width = widths.iter().map(|w| w + 1).sum::<usize>();

    // first header line
    write!(out, "#{:<W$} ", "target", W = target_width)?;
    write!(out, "{:<W$} ", "family", W = fam_width)?;
    write!(out, "{:<W$} ", "%id", W = 3)?;
    tuples
        .iter()
        .zip(&widths)
        .try_for_each(|((_, s, _), w)| write!(out, "{s:<W$} ", W = w))?;
    writeln!(out)?;

    // variable header lines
    for r in 0..n_rows {
        write!(out, "#{:W$} ", "", W = target_width)?;
        write!(out, "{:W$} ", "", W = fam_width)?;
        write!(out, "{:W$} ", "", W = 3)?;
        tuples.iter().zip(&widths).try_for_each(|((p, _, _), w)| {
            let val = p.get(r).map_or("", |v| v);
            write!(out, "{val:<W$} ", W = w)
        })?;
        writeln!(out)?;
    }

    writeln!(
        out,
        "#{}",
        "-".repeat(total_width + target_width + fam_width + 3)
    )?;

    // entries
    for entry in bm.entries {
        write!(out, "{:<W$} ", entry.target, W = target_width)?;
        write!(out, " {:<W$} ", entry.family, W = fam_width)?;
        write!(out, " {:>W$}% ", entry.pid, W = 2)?;
        tuples
            .iter()
            .zip(&widths)
            .try_for_each(|((_, _, m), w)| match m.get(&entry.target) {
                Some(val) => {
                    if args.e_value {
                        write!(out, "{val:<W$.1e} ", W = w)
                    } else {
                        write!(out, "{val:<W$.1} ", W = w)
                    }
                }
                None => write!(out, "{:<W$} ", "-", W = w),
            })?;

        writeln!(out)?;
    }

    Ok(())
}

fn score(args: ScoreArgs) -> anyhow::Result<()> {
    let results = args.which.set().run_dir().join("results");

    let read = |name: &str| -> anyhow::Result<_> {
        let path = manifest::table_path(&results, name, "");
        Ok(bioio::tbl::nail::NailTable::parse(
            File::open(&path).with_context(|| format!("failed to open {}", path.display()))?,
            "",
        )?
        .to_map())
    };

    let full_tbl = read(&args.full)?;
    let sparse_tbl = read(&args.sparse)?;

    let intersection = full_tbl
        .keys()
        .filter(|k| sparse_tbl.contains_key(*k))
        .collect::<Vec<_>>();

    let figures = args.which.out_dir();
    std::fs::create_dir_all(&figures)?;

    let mut out = BufWriter::new(File::create(figures.join("score.txt"))?);
    for k in intersection {
        let f = full_tbl.get(k).unwrap();
        let s = sparse_tbl.get(k).unwrap();
        let x = f.score;
        let y = s.score;
        writeln!(out, "{x:.1},{y:.1}")?;
    }

    Ok(())
}

fn cells(args: CellsArgs) -> anyhow::Result<()> {
    let set = args.which.set();
    let table = manifest::table_path(&set.run_dir().join("results"), &args.run, "");

    let tbl = bioio::tbl::nail::NailTable::parse(
        File::open(&table).with_context(|| format!("failed to open {}", table.display()))?,
        "",
    )?;

    // read out of the files rather than shelled out to hmmstat and
    // esl-seqstat: neither is a dependency this benchmark declares, and both
    // were being found on PATH rather than through `tools`
    let query_lens = bioio::hmm::lengths(set.query_hmm())?;

    let mut target_lens: HashMap<String, usize> = HashMap::new();
    let mut reader = bioio::fasta::Reader::from_path(set.target_fa())?;
    while let Some(rec) = reader.next_record()? {
        target_lens.insert(rec.name.clone(), rec.seq.len());
    }

    let figures = args.which.out_dir();
    std::fs::create_dir_all(&figures)?;

    let mut true_out = BufWriter::new(File::create(figures.join("cells.true.txt"))?);
    let mut decoy_out = BufWriter::new(File::create(figures.join("cells.decoy.txt"))?);

    tbl.hits.iter().try_for_each(|h| -> anyhow::Result<()> {
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

        let x = (qlen * tlen) as f64;
        let y = h.cell_frac;

        if h.target.starts_with("decoy") {
            writeln!(decoy_out, "{x},{y}")?;
        } else if h.query == intended_query {
            writeln!(true_out, "{x},{y}")?;
        }

        Ok(())
    })
}

fn recall(args: RecallArgs) -> anyhow::Result<()> {
    let start = std::time::Instant::now();

    let set = args.which.set();
    let benchmark = Benchmark::new(set.benchmark_tbl())?;
    let data = RecallData::new(&set.run_dir(), &benchmark)?;

    let figures = args.which.out_dir();
    std::fs::create_dir_all(&figures)?;

    let mut roc_path = File::create(figures.join("roc.txt"))?;
    let mut pid_path = File::create(figures.join("pid.txt"))?;
    let mut runtime_path = File::create(figures.join("time.txt"))?;

    data.write_roc(&mut roc_path)?;
    data.write_pid(&mut pid_path)?;
    data.write_runtime(&mut runtime_path)?;

    println!("recall data took: {:?}", start.elapsed());
    Ok(())
}

struct BenchmarkEntry {
    pid: usize,
    target: String,
    query: String,
    family: String,
}

struct Benchmark {
    entries: Vec<BenchmarkEntry>,
    idx_by_pid: HashMap<usize, Vec<usize>>,
    idx_by_target: HashMap<String, usize>,
    idx_by_query: HashMap<String, Vec<usize>>,
    idx_by_family: HashMap<String, Vec<usize>>,
}

impl Benchmark {
    fn new<P: AsRef<Path>>(tbl_path: P) -> anyhow::Result<Self> {
        let tbl_reader = BufReader::new(File::open(tbl_path)?);

        let mut entries = vec![];
        let mut idx_by_pid: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut idx_by_target: HashMap<String, usize> = HashMap::new();
        let mut idx_by_query: HashMap<String, Vec<usize>> = HashMap::new();
        let mut idx_by_family: HashMap<String, Vec<usize>> = HashMap::new();

        for line in tbl_reader
            .lines()
            .map_while(Result::ok)
            .filter(|l| !l.starts_with('#'))
        {
            let tokens: Vec<&str> = line.split_whitespace().collect();

            let pid = tokens[0]
                .strip_suffix('%')
                .context("pid entry missing % symbol")?
                .parse::<usize>()
                .context("")?;
            let family = tokens[1].to_string();
            let target = tokens[2].to_string();
            let query = tokens[3].to_string();

            let entry = BenchmarkEntry {
                pid,
                target: target.clone(),
                query: query.clone(),
                family: family.clone(),
            };

            let at = entries.len();
            idx_by_pid.entry(pid).or_default().push(at);
            idx_by_target.insert(target, at);
            idx_by_query.entry(query).or_default().push(at);
            idx_by_family.entry(family).or_default().push(at);

            entries.push(entry);
        }

        Ok(Self {
            entries,
            idx_by_pid,
            idx_by_target,
            idx_by_query,
            idx_by_family,
        })
    }
}

#[allow(dead_code)]
impl Benchmark {
    fn entries_by_pid(&self, pid: usize) -> Vec<&BenchmarkEntry> {
        self.idx_by_pid
            .get(&pid)
            .unwrap()
            .iter()
            .map(|i| &self.entries[*i])
            .collect()
    }

    fn entry_by_target(&self, target: &str) -> &BenchmarkEntry {
        &self.entries[*self.idx_by_target.get(target).unwrap()]
    }

    fn entries_by_query(&self, query: &str) -> Vec<&BenchmarkEntry> {
        self.idx_by_query
            .get(query)
            .unwrap()
            .iter()
            .map(|i| &self.entries[*i])
            .collect()
    }

    fn entries_by_family(&self, family: &str) -> Vec<&BenchmarkEntry> {
        self.idx_by_family
            .get(family)
            .unwrap()
            .iter()
            .map(|i| &self.entries[*i])
            .collect()
    }
}

#[derive(Clone, Copy, PartialEq)]
enum SearchType {
    Profile,
    Consensus,
    Sequence,
}

impl SearchType {
    /// Off the `mode` field, which the run records rather than the filename
    /// spelling it.
    fn parse(mode: &str) -> anyhow::Result<SearchType> {
        match mode {
            "prf" => Ok(SearchType::Profile),
            "cons" => Ok(SearchType::Consensus),
            "seq" => Ok(SearchType::Sequence),
            other => bail!("unknown search mode {other:?} in manifest.tbl"),
        }
    }
}

/// One finished run: what it was called, how to read its table, which query it
/// searched, and what it cost.
struct Run {
    name: String,
    tool: String,
    mode: SearchType,
    wall_s: f32,
    table: PathBuf,
}

impl Run {
    /// The hits it reported. Which reader to use comes off the `tool` field --
    /// mmseqs, last, blast and diamond all write blast's tabular format,
    /// whatever wrote it.
    fn hits(&self) -> anyhow::Result<HitTable> {
        let file = BufReader::new(
            File::open(&self.table)
                .with_context(|| format!("failed to open {}", self.table.display()))?,
        );

        match self.tool.as_str() {
            "nail" => HitTable::parse::<_, NailTable>(file, &self.name),
            "hmmer" | "phmmer" => HitTable::parse::<_, HmmerTable>(file, &self.name),
            _ => HitTable::parse::<_, BlastTable>(file, &self.name),
        }
        .with_context(|| format!("failed to read {}", self.table.display()))
    }
}

/// The runs a pipeline finished, in the order it declared them.
///
/// Read out of `manifest.tbl` rather than by globbing the results directory,
/// which is what keeps the runs table itself from looking like a hit table and
/// what makes a run's tool and mode facts rather than guesses.
///
/// Several rows can share a name: mmseqs' search and its conversion are one
/// run, and psiblast's per-family calls are one run run a family at a time.
/// [`Wall`] adds those up, and takes the longest rather than the sum of the
/// hmmer parts, which overlap.
fn runs(run_dir: &Path) -> anyhow::Result<Vec<Run>> {
    let manifest = Manifest::read(&run_dir.join("manifest.tbl"))?;
    let results = run_dir.join("results");

    let failed: Vec<&str> = manifest
        .failed()
        .filter_map(|row| row.get(manifest::NAME))
        .collect();
    if !failed.is_empty() {
        eprintln!(
            "warning: leaving out {} command(s) that did not finish: {}",
            failed.len(),
            failed.join(", ")
        );
    }

    let mut out: Vec<Run> = Vec::new();
    let mut walls: Vec<Wall> = Vec::new();
    let mut at: HashMap<String, usize> = HashMap::new();

    for row in manifest.runs() {
        let name = row.get(manifest::NAME).expect("runs() filters on name");

        let i = match at.get(name) {
            Some(&i) => i,
            None => {
                let tool = row
                    .get(manifest::TOOL)
                    .with_context(|| format!("run {name:?} has no tool"))?;
                let mode = row
                    .get(MODE)
                    .with_context(|| format!("run {name:?} has no mode"))?;

                at.insert(name.to_string(), out.len());
                walls.push(Wall::default());
                out.push(Run {
                    name: name.to_string(),
                    tool: tool.to_string(),
                    mode: SearchType::parse(mode)?,
                    wall_s: 0.0,
                    table: manifest::table_path(&results, name, ""),
                });
                out.len() - 1
            }
        };

        walls[i].add("", row);
    }

    for (run, wall) in out.iter_mut().zip(&walls) {
        run.wall_s = wall.total() as f32;
    }

    anyhow::ensure!(
        !out.is_empty(),
        "no finished runs in {}/manifest.tbl",
        run_dir.display()
    );

    Ok(out)
}

impl Display for SearchType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            SearchType::Profile => "prf",
            SearchType::Consensus => "cons",
            SearchType::Sequence => "seq",
        };
        write!(f, "{s}")
    }
}

#[derive(Clone)]
struct Hit2 {
    target: String,
    query: String,
    target_fam: Option<String>,
    query_fam: String,
    e_value: f64,
    pid: Option<usize>,
}

impl Hit2 {
    fn new(hit: &Hit) -> Self {
        let query_tokens: Vec<&str> = hit.query.split('|').collect();
        let query_family = query_tokens[0].to_string();

        let query = if query_tokens.len() > 1 {
            query_tokens[1].to_string()
        } else {
            query_family.clone()
        };

        let (target, target_family, pid) = if hit.target.starts_with("decoy") {
            (hit.target.clone(), None, None)
        } else {
            let target_tokens: Vec<&str> = hit.target.split('|').collect();
            let target_family = target_tokens[0].to_string();
            let target = target_tokens[1].to_string();
            let pid = target_tokens[2]
                .split('%')
                .next()
                .expect("failed to parse pid: no %")
                .parse::<usize>()
                .expect("failed to parse pid");
            (target, Some(target_family), Some(pid))
        };

        Self {
            target,
            query,
            e_value: hit.e_value,
            target_fam: target_family,
            query_fam: query_family,
            pid,
        }
    }
}

impl std::fmt::Display for Hit2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {:.2e}", self.target, self.query, self.e_value)
    }
}

struct HitTable2 {
    name: String,
    positives: Vec<Hit2>,
    /// The decoys this run's own queries reported, in the benchmark's order.
    ///
    /// What the ROC walks: a false positive only counts against a query that
    /// was actually asked, so the list is rebuilt per benchmark entry rather
    /// than taken as everything the run called a decoy.
    adjusted_decoys: Vec<Hit2>,
}

impl HitTable2 {
    fn new(tbl: &HitTable, bm: &Benchmark, search_type: SearchType) -> Self {
        let mut hits_by_target_query_pair: HashMap<(String, String), Vec<Hit2>> = HashMap::new();
        tbl.hits.iter().map(Hit2::new).for_each(|h| {
            hits_by_target_query_pair
                .entry((h.target.clone(), h.query.clone()))
                .or_default()
                .push(h);
        });

        let filtered_hits: Vec<Hit2> = hits_by_target_query_pair
            .into_values()
            .map(|hits| hits.into_iter().min_by(e_value_cmp).expect("empty"))
            .collect();

        let mut positives: Vec<Hit2> = vec![];
        let mut decoys_by_query: HashMap<String, Vec<Hit2>> = HashMap::new();
        filtered_hits.into_iter().for_each(|h| {
            if h.target.starts_with("decoy") {
                decoys_by_query.entry(h.query.clone()).or_default().push(h);
            } else {
                match search_type {
                    SearchType::Profile | SearchType::Consensus => match h.target_fam {
                        Some(ref target_fam) => {
                            if &h.query_fam == target_fam {
                                positives.push(h.clone())
                            }
                        }
                        _ => panic!("hit has no target family"),
                    },
                    SearchType::Sequence => {
                        let paired_query = &bm.entry_by_target(&h.target).query;
                        if &h.query == paired_query {
                            positives.push(h.clone())
                        }
                    }
                }
            }
        });

        let mut adjusted_decoys: Vec<Hit2> = bm
            .entries
            .iter()
            .flat_map(|entry| {
                let name = match search_type {
                    SearchType::Profile => &entry.family,
                    SearchType::Consensus => &format!("{}-consensus", entry.family),
                    SearchType::Sequence => &entry.query,
                };
                match decoys_by_query.get(name) {
                    Some(hits) => hits.clone(),
                    None => vec![],
                }
            })
            .collect();

        positives.sort_by(e_value_cmp);
        adjusted_decoys.sort_by(e_value_cmp);

        Self {
            name: tbl.name.clone(),
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
    fn new(run_dir: &Path, bm: &Benchmark) -> anyhow::Result<Self> {
        let mut pid_bin_tot_cnts = vec![];
        bm.entries.iter().map(|e| e.pid).for_each(|pid| {
            if pid >= pid_bin_tot_cnts.len() {
                pid_bin_tot_cnts.resize(pid + 1, 0);
            }
            pid_bin_tot_cnts[pid] += 1;
        });

        let mut tables = vec![];
        let mut times = vec![];

        for run in runs(run_dir)? {
            tables.push(HitTable2::new(&run.hits()?, bm, run.mode));
            times.push(run.wall_s);
        }

        Ok(Self {
            tables,
            times,
            bin_sizes: pid_bin_tot_cnts,
            positive_cnt: bm.entries.len(),
        })
    }

    fn write_pid<W: Write>(&self, out: &mut W) -> anyhow::Result<()> {
        writeln!(out, "fpr {FIXED_FPR}")?;

        write!(out, "bins")?;
        self.bin_sizes
            .iter()
            .enumerate()
            .try_for_each(|p| write!(out, ",({}, {})", p.0, p.1))?;
        writeln!(out)?;

        let decoy_cnt = (self.positive_cnt as f32 * FIXED_FPR).ceil() as usize;
        self.tables.iter().try_for_each(|tbl| {
            write!(out, "{}", tbl.name)?;

            let e_value_threshold = match tbl.adjusted_decoys.get(decoy_cnt + 1) {
                Some(hit) => hit.e_value,
                None => {
                    println!(
                        "warning: not enough decoys to produce E-value threshold for: {}",
                        tbl.name
                    );
                    f64::INFINITY
                }
            };
            let mut bin_cnts = vec![0usize; self.bin_sizes.len()];
            tbl.positives
                .iter()
                .filter(|h| h.e_value <= e_value_threshold)
                .for_each(|h| {
                    bin_cnts[h.pid.expect("positive hit has no pid")] += 1;
                });

            bin_cnts
                .into_iter()
                .zip(self.bin_sizes.iter())
                .enumerate()
                .try_for_each(|(pid, (cnt, &sz))| {
                    assert!(
                        cnt <= sz,
                        "error: (bin count > bin size): {pid}% | {cnt} > {sz}"
                    );
                    if sz > 0 {
                        write!(
                            out,
                            ",({}, {:.p$})",
                            pid,
                            (cnt as f64 / sz as f64),
                            p = PRECISION
                        )
                    } else {
                        Ok(())
                    }
                })?;
            writeln!(out)
        })?;
        Ok(())
    }

    fn write_roc<W: Write>(&self, out: &mut W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(out);

        let p = 10.0f64.powi(PRECISION as i32);

        self.tables.iter().try_for_each(|tbl| {
            write!(out, "{}", tbl.name)?;

            let mut e_values: Vec<f64> = tbl.adjusted_decoys.iter().map(|h| h.e_value).collect();
            e_values.push(f64::INFINITY);

            let mut counts = vec![];
            let mut last_cnt = 0usize;
            for e in e_values.into_iter() {
                let cnt = tbl.positives[last_cnt..]
                    .iter()
                    .take_while(|h| h.e_value < e)
                    .count();

                last_cnt += cnt;
                counts.push(last_cnt);
            }

            let all_points: Vec<(f64, f64)> = counts
                .into_iter()
                .enumerate()
                .map(|(i, c)| {
                    (
                        i as f64 / self.positive_cnt as f64,
                        c as f64 / self.positive_cnt as f64,
                    )
                })
                .map(|(x, y)| ((x * p).round() / p, (y * p).round() / p))
                .collect();

            let mut points: Vec<(f64, f64)> = vec![all_points[0]];
            all_points
                .windows(2)
                .filter(|p| p[0].1 != p[1].1)
                .for_each(|p| {
                    points.push(p[0]);
                    points.push(p[1]);
                });

            points.dedup();

            points
                .iter()
                .try_for_each(|(x, y)| write!(out, ",({x:.p$}, {y:.p$})", p = PRECISION))?;

            writeln!(out)
        })?;
        Ok(())
    }

    fn write_runtime<W: Write>(&self, out: &mut W) -> anyhow::Result<()> {
        writeln!(out, "fpr {FIXED_FPR}")?;

        let decoy_cnt = (self.positive_cnt as f32 * FIXED_FPR).ceil() as usize;
        self.tables
            .iter()
            .zip(self.times.iter())
            .try_for_each(|(tbl, time)| {
                let e_value_threshold = match tbl.adjusted_decoys.get(decoy_cnt + 1) {
                    Some(hit) => hit.e_value,
                    None => {
                        println!(
                            "warning: not enough decoys to produce E-value threshold for: {}",
                            tbl.name
                        );
                        f64::INFINITY
                    }
                };

                let count = tbl
                    .positives
                    .iter()
                    .take_while(|h| h.e_value < e_value_threshold)
                    .count();

                let recall = count as f64 / self.positive_cnt as f64;

                writeln!(out, "{},({:.4},{:.4})", tbl.name, time, recall)
            })?;
        Ok(())
    }
}
