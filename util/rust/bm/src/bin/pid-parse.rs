use std::{
    collections::HashMap,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use bioio::tbl::{BlastTable, Hit, HitTable, HmmerTable, NailTable};

use anyhow::Context;
use clap::{Parser, Subcommand};
use glob::glob;
use regex::Regex;

const PRECISION: usize = 4;
const FIXED_FPR: f32 = 0.01;

trait Float: PartialOrd {}
impl Float for f32 {}
impl Float for f64 {}

fn float_cmp<F: Float>(a: &F, b: &F) -> std::cmp::Ordering {
    a.partial_cmp(b).expect("NaN encountered in float cmp")
}

fn e_value_cmp(a: &Hit2, b: &Hit2) -> std::cmp::Ordering {
    a.e_value
        .partial_cmp(&b.e_value)
        .expect("NaN encountered in E-value cmp")
}

#[derive(Parser)]
struct RecallArgs {
    #[arg(value_name = "benchmark.tbl")]
    benchmark_tbl: PathBuf,

    #[arg(value_name = "results/")]
    results_dir: PathBuf,

    #[arg(value_name = "dir", default_value = "./figures")]
    out_dir: PathBuf,
}

#[derive(Parser)]
struct CellsArgs {
    #[arg(value_name = "nail.tbl")]
    nail_tbl: PathBuf,

    #[arg(value_name = "query.hmm")]
    query: PathBuf,

    #[arg(value_name = "target.fa")]
    target: PathBuf,

    #[arg(value_name = "dir", default_value = "./figures")]
    out_dir: PathBuf,
}

#[derive(Parser)]
struct ScoreArgs {
    #[arg(value_name = "nail.full.tbl")]
    full_tbl: PathBuf,

    #[arg(value_name = "nail.sparse.tbl")]
    sparse_tbl: PathBuf,

    #[arg(value_name = "dir", default_value = "./figures")]
    out_dir: PathBuf,
}

#[derive(Parser)]
struct TableArgs {
    #[arg(short, long, value_name = "benchmark.tbl")]
    benchmark_tbl: PathBuf,

    #[arg(short, long, value_name = "results/")]
    results_dir: PathBuf,

    #[arg(short, long, value_name = "out.tbl", default_value = "results.tbl")]
    out_path: PathBuf,

    #[arg(short, long, default_value_t = false)]
    e_value: bool,

    #[arg(
        short,
        long,
        default_value_t = 6,
        default_value_if("e_value", "true", "9")
    )]
    min_width: usize,
}

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    Recall(RecallArgs),
    Cells(CellsArgs),
    Score(ScoreArgs),
    Table(TableArgs),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().cmd {
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
    let bm = Benchmark::new(args.benchmark_tbl).context("failed to open benchmark.tbl")?;

    let mut out = BufWriter::new(File::create(args.out_path)?);

    let mut tuples = vec![];

    for path in glob(
        args.results_dir
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    {
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("failed to produce file stem")?;

        let (prefix, search_type) = stem.rsplit_once(".").context("failed to split prefix")?;
        let prefix_tokens = prefix.split("-").map(|s| s.to_string()).collect::<Vec<_>>();

        let file = BufReader::new(File::open(&path)?);

        let tbl = match prefix_tokens.first().context("no prefix tokens")?.as_str() {
            "nail" => HitTable::parse::<_, NailTable>(file, prefix)?,
            "hmmer" => HitTable::parse::<_, HmmerTable>(file, prefix)?,
            "mmseqs" => HitTable::parse::<_, BlastTable>(file, prefix)?,
            "last" => HitTable::parse::<_, BlastTable>(file, prefix)?,
            "blast" => HitTable::parse::<_, BlastTable>(file, prefix)?,
            "diamond" => HitTable::parse::<_, BlastTable>(file, prefix)?,
            _ => continue,
        };

        fn true_hit_filter(hit: &Hit) -> bool {
            if hit.target.starts_with("decoy") {
                return false;
            }

            let q = hit.query.split('|').next().unwrap();
            let t = hit.target.split('|').next().unwrap();

            q == t
        }

        let hits = tbl
            .hits
            .into_iter()
            .filter(true_hit_filter)
            .collect::<Vec<_>>();

        let map: HashMap<String, f64> = if search_type == "prf" {
            hits.into_iter()
                .map(|hit| {
                    let t = hit.target.split('|').nth(1).unwrap().to_string();

                    if args.e_value {
                        (t, hit.e_value)
                    } else {
                        (t, hit.score as f64)
                    }
                })
                .collect()
        } else {
            // if we have seq queries, retain the single highest score
            let mut best: HashMap<String, f64> = HashMap::new();
            for hit in hits {
                let t = hit.target.split('|').nth(1).unwrap().to_string();

                if args.e_value {
                    best.entry(t)
                        .and_modify(|score| *score = score.min(hit.e_value))
                        .or_insert(hit.e_value);
                } else {
                    best.entry(t)
                        .and_modify(|score| *score = score.max(hit.score as f64))
                        .or_insert(hit.score as f64);
                }
            }
            best
        };

        tuples.push((prefix_tokens, search_type.to_string(), map));
    }

    // ---

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

    // ---

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
    let full_tbl = bioio::tbl::nail::NailTable::parse(
        File::open(&args.full_tbl)
            .with_context(|| format!("failed to open: {:?}", args.full_tbl))?,
        "",
    )?
    .to_map();

    let sparse_tbl = bioio::tbl::nail::NailTable::parse(
        File::open(&args.sparse_tbl)
            .with_context(|| format!("failed to open: {:?}", args.sparse_tbl))?,
        "",
    )?
    .to_map();

    let intersection = full_tbl
        .keys()
        .filter(|k| sparse_tbl.contains_key(*k))
        .collect::<Vec<_>>();

    let figures = args.out_dir;
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
    let tbl = bioio::tbl::nail::NailTable::parse(File::open(args.nail_tbl)?, "")?;

    let query_lens: HashMap<String, usize> = Command::new("hmmstat")
        .arg(args.query)
        .output()?
        .stdout
        .lines()
        .filter_map(|line| {
            let line = line.ok()?;
            if line.starts_with('#') {
                None
            } else {
                let tokens: Vec<&str> = line.split_whitespace().collect();
                Some((tokens[1].to_string(), tokens[5].parse::<usize>().ok()?))
            }
        })
        .collect();

    let target_lens: HashMap<String, usize> = Command::new("esl-seqstat")
        .arg("-a")
        .arg(args.target)
        .output()?
        .stdout
        .lines()
        .filter_map(|line| {
            let line = line.ok()?;
            if line.starts_with('=') {
                let tokens: Vec<&str> = line.split_whitespace().collect();
                Some((tokens[1].to_string(), tokens[2].parse::<usize>().ok()?))
            } else {
                None
            }
        })
        .collect();

    let figures = args.out_dir;
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

    let benchmark = Benchmark::new(args.benchmark_tbl).context("failed to open benchmark.tbl")?;

    let data = RecallData::new(args.results_dir, &benchmark).context("failed to get data")?;

    let figures = args.out_dir;
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
    target_cnt: usize,
    entries: Vec<BenchmarkEntry>,
    idx_by_pid: HashMap<usize, Vec<usize>>,
    idx_by_target: HashMap<String, usize>,
    idx_by_query: HashMap<String, Vec<usize>>,
    idx_by_family: HashMap<String, Vec<usize>>,
}

impl Benchmark {
    fn new<P: AsRef<Path>>(tbl_path: P) -> anyhow::Result<Self> {
        let tbl_reader = BufReader::new(File::open(tbl_path)?);

        let mut target_cnt = 0;
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

            idx_by_pid.entry(pid).or_default().push(target_cnt);
            idx_by_target.insert(target, target_cnt);
            idx_by_query.entry(query).or_default().push(target_cnt);
            idx_by_family.entry(family).or_default().push(target_cnt);

            entries.push(entry);
            target_cnt += 1;
        }

        Ok(Self {
            target_cnt,
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

#[derive(PartialEq)]
enum SearchType {
    Profile,
    Consensus,
    Sequence,
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
    search_type: SearchType,
    positives: Vec<Hit2>,
    decoys: Vec<Hit2>,
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
        let mut decoys: Vec<Hit2> = vec![];
        let mut decoys_by_query: HashMap<String, Vec<Hit2>> = HashMap::new();
        filtered_hits.into_iter().for_each(|h| {
            if h.target.starts_with("decoy") {
                let query_decoys = decoys_by_query.entry(h.query.clone()).or_default();
                query_decoys.push(h.clone());
                decoys.push(h.clone());
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
        decoys.sort_by(e_value_cmp);
        adjusted_decoys.sort_by(e_value_cmp);

        Self {
            name: tbl.name.clone(),
            search_type,
            positives,
            decoys,
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
    fn new<P: AsRef<Path>>(results_dir: P, bm: &Benchmark) -> anyhow::Result<Self> {
        let results_dir = PathBuf::from(results_dir.as_ref());
        let mut times = vec![];
        let mut tables = vec![];

        let mut pid_bin_tot_cnts = vec![];
        bm.entries.iter().map(|e| e.pid).for_each(|pid| {
            if pid >= pid_bin_tot_cnts.len() {
                pid_bin_tot_cnts.resize(pid + 1, 0);
            }
            pid_bin_tot_cnts[pid] += 1;
        });

        let time_pattern = Regex::new(r"Elapsed.*\): (?P<time>.*)$").unwrap();

        for path in glob(
            results_dir
                .join("*.tbl")
                .to_str()
                .context("invalid *.tbl glob")?,
        )?
        .filter_map(Result::ok)
        {
            let file = BufReader::new(File::open(&path)?);

            let stem_tokens: Vec<&str> = path
                .file_stem()
                .and_then(|s| s.to_str())
                .context("invalid path")?
                .split('.')
                .collect();

            assert!(stem_tokens.len() >= 2);

            let prefix = stem_tokens[..stem_tokens.len() - 1].join(".");
            let search_type = match *stem_tokens.last().unwrap() {
                "cons" => SearchType::Consensus,
                "prf" => SearchType::Profile,
                "seq" => SearchType::Sequence,
                _ => panic!("unknown search type"),
            };

            let name = format!("{prefix}.{search_type}");

            let tbl = match prefix {
                s if s.starts_with("hmmer") => HitTable::parse::<_, HmmerTable>(file, &name),
                s if s.starts_with("nail") => HitTable::parse::<_, NailTable>(file, &name),
                _ => HitTable::parse::<_, BlastTable>(file, &name),
            }
            .map(|tbl| HitTable2::new(&tbl, bm, search_type))?;
            tables.push(tbl);

            // the time file will have the same prefix as the tbl
            let time_path = path.with_extension("time");
            let time: f32 = std::fs::read_to_string(time_path)?
                .lines()
                .filter_map(|l| time_pattern.captures(l))
                .map(|c| {
                    c.name("time")
                        .unwrap()
                        .as_str()
                        .split(':')
                        .map(|t| t.parse::<f32>().unwrap())
                        .rev()
                        .enumerate()
                        .map(|(i, t)| t * 60.0_f32.powf(i as f32))
                        .sum()
                })
                .max_by(float_cmp)
                .expect("no times found");

            times.push(time);
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
            // println!("{} {:.3e}", tbl.name, e_value_threshold);
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
