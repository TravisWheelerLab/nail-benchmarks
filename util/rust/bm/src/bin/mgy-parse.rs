use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use bioio::tbl::{hmmer::HmmerDomainTable, BlastTable, HitTable, HmmerTable, NailTable};

use clap::{Parser, Subcommand};
use glob::glob;
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    ThreadPoolBuilder,
};
use regex::Regex;

trait Float: PartialOrd {}
impl Float for f32 {}
impl Float for f64 {}

fn float_cmp<F: Float>(a: &F, b: &F) -> std::cmp::Ordering {
    a.partial_cmp(b).expect("NaN encountered in float cmp")
}

fn hit_cmp(a: &bioio::tbl::Hit, b: &bioio::tbl::Hit) -> std::cmp::Ordering {
    float_cmp(&a.e_value, &b.e_value)
}

fn set_threads(num_threads: usize) -> anyhow::Result<()> {
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("failed to build rayon global threadpool")
}

fn p_value(score: f64, lambda: f64, tau: f64) -> f64 {
    (-lambda * (score - tau)).exp()
}

struct HmmGumbel {
    ga_sc: f64,
    ga_p: f64,
    tau: f64,
    lambda: f64,
}

impl HmmGumbel {
    fn p_value(&self, score: f64) -> f64 {
        (-self.lambda * (score - self.tau)).exp()
    }
}

fn target_db_size(target_path: impl AsRef<Path>) -> anyhow::Result<f64> {
    let reader = BufReader::new(File::open(target_path.as_ref())?);
    let mut z = 0.0;

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            z += 1.0;
        }
    }

    Ok(z)
}

fn hmms(hmm_path: impl AsRef<Path>) -> anyhow::Result<HashMap<String, HmmGumbel>> {
    let reader = BufReader::new(File::open(hmm_path.as_ref())?);

    let mut names = vec![];
    let mut gathering_thresholds = vec![];
    let mut gumbels = vec![];

    for line in reader.lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix("NAME") {
            names.push(rest.split_whitespace().collect::<String>());
        }

        if let Some(rest) = line.strip_prefix("GA") {
            let x = rest
                .split_whitespace()
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>();

            gathering_thresholds.push((x[0], x[1]));
        }

        if let Some(rest) = line.strip_prefix("STATS LOCAL FORWARD") {
            let x = rest
                .split_whitespace()
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>();

            gumbels.push((x[0], x[1]))
        }
    }

    assert_eq!(names.len(), gathering_thresholds.len());
    assert_eq!(names.len(), gumbels.len());

    Ok(names
        .into_iter()
        .enumerate()
        .map(|(i, n)| {
            let ga_sc = gathering_thresholds[i].0;
            let tau = gumbels[i].0;
            let lambda = gumbels[i].1;
            let ga_p = (-lambda * (ga_sc - tau)).exp();
            (
                n,
                HmmGumbel {
                    ga_sc,
                    ga_p,
                    tau,
                    lambda,
                },
            )
        })
        .collect())
}

#[derive(Default)]
pub struct Record {
    query: String,
    target: String,
    ga_score: f32,
    ga_p_value: f64,
    nail_cloud_score: Option<f32>,
    nail_cloud_p_value: Option<f64>,
    nail_score: Option<f32>,
    nail_p_value: Option<f64>,
    mmseqs_score: Option<f32>,
    mmseqs_p_value: Option<f64>,
    hmmer_score: Option<f32>,
    hmmer_p_value: Option<f64>,
    dom_score_sum: Option<f32>,
    dom_score_max: Option<f32>,
    dom_sig_cnt: Option<usize>,
    dom_scores: Option<Vec<f32>>,
}

fn score_fmt(f: Option<f32>, w: usize) -> String {
    match f {
        Some(s) => format!("{s:^W$.1}", W = w),
        None => format!("{:^W$}", "-", W = w),
    }
}

fn p_value_fmt(f: Option<f64>, w: usize) -> String {
    match f {
        Some(s) => format!("{s:^W$.2e}", W = w),
        None => format!("{:^W$}", "-", W = w),
    }
}

fn int_fmt(f: Option<usize>, w: usize) -> String {
    match f {
        Some(s) => format!("{s:^W$}", W = w),
        None => format!("{:^W$}", "-", W = w),
    }
}

const HEADER: [&str; 3] = [
    "| GA |    GA     |  nail  |    nail   |  nail  |   nail    | mmseqs |  mmseqs   | hmmer |  hmmer    | dom | dom | sig |",
    "| sc |  P-value  | cld sc |cld P-value|   sc   |  P-value  |   sc   |  P-value  |  sc   |  P-value  | sum | max | dom | dom scores",
    " ---- ----------- -------- ----------- -------- ----------- -------- ----------- ------- ----------- ----- ----- ----- ------------",
];

fn write_header(out: &mut impl Write, recs: &[Record]) -> anyhow::Result<()> {
    let q_max = recs.iter().map(|r| r.query.len()).max().unwrap();
    let t_max = recs.iter().map(|r| r.target.len()).max().unwrap();

    writeln!(
        out,
        "#{:W1$} {:W2$}{}",
        " ",
        " ",
        HEADER[0],
        W1 = q_max - 1,
        W2 = t_max
    )?;
    writeln!(
        out,
        "#{:^W1$} {:^W2$}{}",
        "query",
        "target",
        HEADER[1],
        W1 = q_max - 1,
        W2 = t_max,
    )?;
    writeln!(
        out,
        "#{:W1$} {:W2$}{}",
        "-".repeat(q_max - 1),
        "-".repeat(t_max),
        HEADER[2],
        W1 = q_max - 1,
        W2 = t_max,
    )?;

    Ok(())
}

fn write_records(out: &mut impl Write, recs: &[Record], w: &[usize]) -> anyhow::Result<()> {
    let q_max = recs.iter().map(|r| r.query.len()).max().unwrap();
    let t_max = recs.iter().map(|r| r.target.len()).max().unwrap();

    for r in recs.iter() {
        writeln!(
            out,
            "{:W1$} {:W2$} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            r.query,
            r.target,
            score_fmt(Some(r.ga_score), w[0]),
            p_value_fmt(Some(r.ga_p_value), w[1]),
            score_fmt(r.nail_cloud_score, w[2]),
            p_value_fmt(r.nail_cloud_p_value, w[3]),
            score_fmt(r.nail_score, w[4]),
            p_value_fmt(r.nail_p_value, w[5]),
            score_fmt(r.mmseqs_score, w[6]),
            p_value_fmt(r.mmseqs_p_value, w[7]),
            score_fmt(r.hmmer_score, w[8]),
            p_value_fmt(r.hmmer_p_value, w[9]),
            score_fmt(r.dom_score_sum, w[10]),
            score_fmt(r.dom_score_max, w[11]),
            int_fmt(r.dom_sig_cnt, w[12]),
            match r.dom_scores {
                Some(ref v) => v
                    .iter()
                    .map(|f| format!("{f:.1}"))
                    .collect::<Vec<_>>()
                    .join(","),
                None => "-".to_string(),
            },
            W1 = q_max,
            W2 = t_max
        )?;
    }
    Ok(())
}

#[derive(Subcommand)]
enum SubCommands {
    Standard(StandardArgs),
    Rev(RevArgs),
    Params(ParamsArgs),
}

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    pub command: SubCommands,
}

#[derive(Parser)]
struct StandardArgs {
    #[arg(value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(value_name = "target.fa")]
    target_path: PathBuf,

    #[arg(value_name = "hmmer.tbl")]
    hmmer_tbl_path: PathBuf,

    #[arg(value_name = "hmmer.domtbl")]
    hmmer_domtbl_path: PathBuf,

    #[arg(value_name = "nail.tbl")]
    nail_tbl_path: PathBuf,

    #[arg(long, value_name = "nail.stats")]
    nail_stats_path: Option<PathBuf>,

    #[arg(value_name = "mmseqs.tbl")]
    mmseqs_tbl_path: PathBuf,

    #[arg(value_name = "out.tbl")]
    out_tbl_path: PathBuf,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,
}

#[derive(Parser)]
struct ParamsArgs {
    #[arg(value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(value_name = "target.fa")]
    target_path: PathBuf,

    #[arg(value_name = "hmmer.tbl")]
    hmmer_tbl_path: PathBuf,

    #[arg(value_name = "nail/")]
    nail_dir: PathBuf,

    #[arg(value_name = "mmseqs/")]
    mmseqs_dir: PathBuf,

    #[arg(value_name = "out.txt")]
    out_path: PathBuf,
}

#[derive(Parser)]
struct RevArgs {
    #[arg(value_name = "nail.tbl")]
    nail_dir: PathBuf,

    #[arg(value_name = "mmseqs.tbl")]
    mmseqs_dir: PathBuf,

    #[arg(value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(value_name = "out/")]
    out_dir: PathBuf,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        SubCommands::Standard(args) => standard(args),
        SubCommands::Params(args) => params(args),
        SubCommands::Rev(args) => rev(args),
    }
}

fn rev(args: RevArgs) -> anyhow::Result<()> {
    let start = std::time::Instant::now();

    let figures = args.out_dir;
    std::fs::create_dir_all(&figures)?;

    set_threads(args.num_threads)?;

    // ---

    let hmms =
        hmms(&args.query_path).with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    // ---

    let mut pieces = glob(
        args.nail_dir
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    .filter_map(|p| {
        p.file_name()
            .and_then(|n| n.to_str())
            .and_then(|s| s.split('.').nth(1))
            .and_then(|i| i.parse::<usize>().ok())
    })
    .collect::<Vec<_>>();

    pieces.sort();
    pieces.dedup();

    // ---

    type DecoyMap = HashMap<String, Vec<bioio::tbl::Hit>>;

    struct RevTables {
        nail: bioio::tbl::HitTable,
        nail_rev: bioio::tbl::HitTable,
        mmseqs: bioio::tbl::HitTable,
        mmseqs_rev: bioio::tbl::HitTable,
    }

    impl RevTables {
        fn new(nail_dir: impl AsRef<Path>, mmseqs_dir: impl AsRef<Path>, idx: usize) -> Self {
            let nail_dir = nail_dir.as_ref();
            let mmseqs_dir = mmseqs_dir.as_ref();

            let path = nail_dir.join(format!("nail.{idx}.prf.tbl"));
            let nail = HitTable::from_path::<_, NailTable>(&path).expect("failed to open {path:?}");

            let path = nail_dir.join(format!("nail.{idx}.rev.prf.tbl"));
            let nail_rev =
                HitTable::from_path::<_, NailTable>(&path).expect("failed to open {path:?}");

            let path = mmseqs_dir.join(format!("mmseqs.{idx}.prf.tbl"));
            let mmseqs =
                HitTable::from_path::<_, BlastTable>(&path).expect("failed to open {path:?}");

            let path = mmseqs_dir.join(format!("mmseqs.{idx}.rev.prf.tbl"));
            let mmseqs_rev =
                HitTable::from_path::<_, BlastTable>(&path).expect("failed to open {path:?}");

            Self {
                nail,
                nail_rev,
                mmseqs,
                mmseqs_rev,
            }
        }
    }

    // type ScoreMap = HashMap<String, f64>;
    // #[derive(Default)]
    // struct DecoyScoreDistData {
    //     nail_min: ScoreMap,
    //     nail_max: ScoreMap,
    //     mmseqs_min: ScoreMap,
    //     mmseqs_max: ScoreMap,
    // }

    // impl DecoyScoreDistData {
    //     fn update(&mut self, tables: RevTables) {
    //         let mut nail_map = tables.nail.to_map();
    //         let mut nail_rev_map = tables.nail_rev.to_map();
    //         let mut mmseqs_map = tables.mmseqs.to_map();
    //         let mut mmseqs_rev_map = tables.mmseqs_rev.to_map();

    //         nail_map.values().for_each(|h| {
    //             self.nail_min
    //                 .entry(h.query.clone())
    //                 .and_modify(|s| *s = s.min(h.score))
    //                 .or_insert(h.score);

    //             self.nail_max
    //                 .entry(h.query.clone())
    //                 .and_modify(|s| *s = s.max(h.score))
    //                 .or_insert(h.score);
    //         });

    //         mmseqs_map.values().for_each(|h| {
    //             self.mmseqs_min
    //                 .entry(h.query.clone())
    //                 .and_modify(|s| *s = s.min(h.score))
    //                 .or_insert(h.score);

    //             self.mmseqs_max
    //                 .entry(h.query.clone())
    //                 .and_modify(|s| *s = s.max(h.score))
    //                 .or_insert(h.score);
    //         });
    //     }

    //     fn merge(&mut self, other: Self) {
    //         for (k, v) in other.nail_min {
    //             self.nail_min
    //                 .entry(k)
    //                 .and_modify(|min| *min = min.min(v))
    //                 .or_insert(v);
    //         }

    //         for (k, v) in other.nail_max {
    //             self.nail_max
    //                 .entry(k)
    //                 .and_modify(|max| *max = max.max(v))
    //                 .or_insert(v);
    //         }

    //         for (k, v) in other.mmseqs_min {
    //             self.mmseqs_min
    //                 .entry(k)
    //                 .and_modify(|min| *min = min.min(v))
    //                 .or_insert(v);
    //         }

    //         for (k, v) in other.mmseqs_max {
    //             self.mmseqs_max
    //                 .entry(k)
    //                 .and_modify(|max| *max = max.max(v))
    //                 .or_insert(v);
    //         }
    //     }
    // }

    #[derive(Default)]
    struct DecoyCutoffData {
        nail_decoys: DecoyMap,
        mmseqs_decoys: DecoyMap,
    }

    impl DecoyCutoffData {
        fn update(&mut self, mut tables: RevTables) {
            // ---
            // filter the real hits by E-value

            tables.nail.hits.retain(|h| h.e_value <= E);
            tables.mmseqs.hits.retain(|h| h.e_value <= E);

            // ---
            // convert to (query, target)-keyed maps for easier comparison

            let nail_map = tables.nail.to_map();
            let mut nail_rev_map = tables.nail_rev.to_map();
            let mmseqs_map = tables.mmseqs.to_map();
            let mut mmseqs_rev_map = tables.mmseqs_rev.to_map();

            // ---
            // retain only reverse hits for pairs that don't
            // remain in the real hits after filtering

            nail_rev_map.retain(|k, _| !nail_map.contains_key(k));
            mmseqs_rev_map.retain(|k, _| !mmseqs_map.contains_key(k));

            // ---

            nail_rev_map
                .into_values()
                .for_each(|h| self.nail_decoys.entry(h.query.clone()).or_default().push(h));

            mmseqs_rev_map.into_values().for_each(|h| {
                self.mmseqs_decoys
                    .entry(h.query.clone())
                    .or_default()
                    .push(h)
            });
        }

        fn merge(&mut self, other: Self) {
            for (k, mut v) in other.nail_decoys {
                self.nail_decoys.entry(k).or_default().append(&mut v);
            }

            for (k, mut v) in other.mmseqs_decoys {
                self.mmseqs_decoys.entry(k).or_default().append(&mut v);
            }
        }
    }

    const E: f64 = 1e-3;
    let mut data = pieces
        .par_iter()
        .fold(DecoyCutoffData::default, |mut data, &idx| {
            let tables = RevTables::new(&args.nail_dir, &args.mmseqs_dir, idx);
            data.update(tables);
            data
        })
        .reduce(DecoyCutoffData::default, |mut d1, d2| {
            d1.merge(d2);
            d1
        });

    // ---

    queries.iter().for_each(|q| {
        data.nail_decoys.entry(q.to_string()).or_default();
        data.mmseqs_decoys.entry(q.to_string()).or_default();
    });

    data.nail_decoys
        .values_mut()
        .for_each(|v| v.sort_by(hit_cmp));
    data.mmseqs_decoys
        .values_mut()
        .for_each(|v| v.sort_by(hit_cmp));

    // ---

    let mut out_cutoffs = BufWriter::new(File::create(figures.join("cutoffs.txt"))?);

    let mut out_dist = BufWriter::new(File::create(figures.join("dist.txt"))?);
    let mut out_ga = BufWriter::new(File::create(figures.join("ga.txt"))?);
    let mut out_cnt = BufWriter::new(File::create(figures.join("count.txt"))?);

    const N_CUTOFF: usize = 5;
    for (q, nail_hits) in data.nail_decoys.iter() {
        let hmm = hmms.get(q).unwrap();
        let mmseqs_hits = data.mmseqs_decoys.get(q).unwrap();

        let nail_scores = nail_hits
            .iter()
            .map(|t| t.score)
            .chain(std::iter::repeat(0.0))
            .take(N_CUTOFF)
            .collect::<Vec<_>>();

        let mmseqs_scores = mmseqs_hits
            .iter()
            .map(|t| t.score)
            .chain(std::iter::repeat(0.0))
            .take(N_CUTOFF)
            .collect::<Vec<_>>();

        writeln!(
            out_cutoffs,
            "{q},(nail,{}),(mmseqs,{})",
            nail_scores
                .iter()
                .map(|s| format!("{s:.1}"))
                .collect::<Vec<_>>()
                .join(","),
            mmseqs_scores
                .iter()
                .map(|s| format!("{s:.1}"))
                .collect::<Vec<_>>()
                .join(","),
        )?;

        // ---

        writeln!(out_dist, "{q},({},{})", nail_hits.len(), mmseqs_hits.len())?;

        writeln!(out_ga, "{:.1},{:.1}", nail_scores[0], hmm.ga_sc)?;

        let diff = hmm.ga_sc - nail_scores[0];
        writeln!(out_cnt, "{},{diff:.1}", nail_hits.len())?;
    }

    // ---

    println!("took {:?}", start.elapsed());
    Ok(())
}

fn params(args: ParamsArgs) -> anyhow::Result<()> {
    let start = std::time::Instant::now();

    let z = target_db_size(args.target_path)?;

    let hmms = hmms(args.query_path)?;

    let hmmer_tbl = HitTable::from_path::<_, HmmerTable>(args.hmmer_tbl_path)?;

    let mut ga_hmmer: HashSet<(&str, &str)> = HashSet::new();

    for h in &hmmer_tbl.hits {
        let hmm = hmms.get(&h.query).unwrap();

        if h.score >= hmm.ga_sc {
            ga_hmmer.insert((&h.query, &h.target));
        }
    }

    let hmmer_cnt = ga_hmmer.len() as f32;

    let time_pattern = Regex::new(r"Elapsed.*\): (?P<time>.*)$").unwrap();

    let mut out = BufWriter::new(File::create(&args.out_path)?);

    for path in glob(
        args.nail_dir
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    {
        let mut cnt = 0.0;

        let tbl = HitTable::from_path::<_, NailTable>(&path)?;
        for h in &tbl.hits {
            let hmm = hmms.get(&h.query).unwrap();

            if ga_hmmer.contains(&(&h.query, &h.target)) && h.score >= hmm.ga_sc {
                cnt += 1.0;
            }
        }

        let frac = cnt / hmmer_cnt;

        let stem_tokens: Vec<&str> = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("invalid path")?
            .split('.')
            .collect();

        let prefix = stem_tokens[..stem_tokens.len() - 1].join(".");
        let search_type = stem_tokens.last().unwrap();

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

        writeln!(out, "{prefix},{search_type},({time:.4},{frac:.4})")?;
    }

    for path in glob(
        args.mmseqs_dir
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    {
        let mut cnt = 0.0;

        let tbl = HitTable::from_path::<_, BlastTable>(&path)?;
        for h in &tbl.hits {
            let hmm = hmms.get(&h.query).unwrap();

            if ga_hmmer.contains(&(&h.query, &h.target)) && h.e_value / z <= hmm.ga_p {
                cnt += 1.0;
            }
        }

        let frac = cnt / hmmer_cnt;

        let stem_tokens: Vec<&str> = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("invalid path")?
            .split('.')
            .collect();

        let prefix = stem_tokens[..stem_tokens.len() - 1].join(".");
        let search_type = stem_tokens.last().unwrap();

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

        writeln!(out, "{prefix},{search_type},({time:.4},{frac:.4})")?;
    }

    println!(
        "{} took {:.2}s",
        args.out_path.to_string_lossy(),
        start.elapsed().as_secs_f32()
    );

    Ok(())
}

fn standard(args: StandardArgs) -> anyhow::Result<()> {
    let start = std::time::Instant::now();

    set_threads(args.num_threads)?;

    let z = target_db_size(args.target_path)?;

    let hmms = hmms(args.query_path)?;

    let nail_tbl = HitTable::from_path::<_, NailTable>(args.nail_tbl_path)?;
    let hmmer_tbl = HitTable::from_path::<_, HmmerTable>(args.hmmer_tbl_path)?;
    let mmseqs_tbl = HitTable::from_path::<_, BlastTable>(args.mmseqs_tbl_path)?;

    let mut ga_any: HashSet<(&str, &str)> = HashSet::new();

    // ---

    for h in &nail_tbl.hits {
        let hmm = hmms.get(&h.query).unwrap();

        if h.score >= hmm.ga_sc {
            ga_any.insert((&h.query, &h.target));
        }
    }

    for h in &hmmer_tbl.hits {
        let hmm = hmms.get(&h.query).unwrap();

        if h.score >= hmm.ga_sc {
            ga_any.insert((&h.query, &h.target));
        }
    }

    for h in &mmseqs_tbl.hits {
        let hmm = hmms.get(&h.query).unwrap();

        if h.e_value / z <= hmm.ga_p {
            ga_any.insert((&h.query, &h.target));
        }
    }

    // ---

    let mut records_by_pair = ga_any
        .iter()
        .map(|(q, t)| {
            let hmm = hmms.get(*q).unwrap();

            (
                (q.to_string(), t.to_string()),
                Record {
                    query: q.to_string(),
                    target: t.to_string(),
                    ga_score: hmm.ga_sc as f32,
                    ga_p_value: hmm.ga_p,
                    ..Default::default()
                },
            )
        })
        .collect::<HashMap<(String, String), Record>>();

    let test: HashSet<_> = ga_any
        .iter()
        .map(|(q, t)| (q.to_string(), t.to_string()))
        .collect();

    let test = std::sync::Arc::new(test);

    // ---

    for h in nail_tbl.hits {
        if let Some(r) = records_by_pair.get_mut(&(h.query, h.target)) {
            r.nail_score = Some(h.score as f32);
            r.nail_p_value = Some(h.e_value / z);
        }
    }

    for h in hmmer_tbl.hits {
        if let Some(r) = records_by_pair.get_mut(&(h.query, h.target)) {
            r.hmmer_score = Some(h.score as f32);
            r.hmmer_p_value = Some(h.e_value / z);
        }
    }

    for h in mmseqs_tbl.hits {
        if let Some(r) = records_by_pair.get_mut(&(h.query, h.target)) {
            r.mmseqs_score = Some(h.score as f32);
            r.mmseqs_p_value = Some(h.e_value / z);
        }
    }

    // --

    if let Some(path) = args.nail_stats_path {
        let reader = BufReader::new(File::open(path)?);

        let mut it = reader.lines();

        while let Ok(batch) = it.by_ref().take(100_000).collect::<Result<Vec<_>, _>>() {
            if batch.is_empty() {
                break;
            }

            let filtered = batch
                .par_iter()
                .filter_map(|line| {
                    let mut it = line.split_whitespace();

                    let q = it.next()?.to_string();
                    let t = it.next()?.to_string();

                    // skip 4 fields
                    it.nth(3)?;

                    let sc = it.next()?.parse::<f32>().ok()?;
                    let p = it.next()?.parse::<f64>().ok()?;

                    let key = (q, t);
                    test.contains(&key).then_some((key, sc, p))
                })
                .collect::<Vec<_>>();

            for (k, sc, p) in filtered {
                if let Some(r) = records_by_pair.get_mut(&k) {
                    r.nail_cloud_score = Some(sc);
                    r.nail_cloud_p_value = Some(p);
                }
            }
        }
    }
    // --

    let dom_tbl = HmmerDomainTable::from_path(args.hmmer_domtbl_path)?;

    for (k, v) in records_by_pair.iter_mut() {
        let hit = match dom_tbl.hits.get(k) {
            Some(d) => d,
            None => continue,
        };

        let mut dom = hit.domains.iter().map(|d| d.score).collect::<Vec<_>>();

        let dom_score_sum = dom.iter().sum::<f32>();
        dom.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let dom_score_max = *dom.last().unwrap();

        dom.reverse();

        let dom_pct = dom.iter().map(|s| s / dom_score_max).collect::<Vec<_>>();
        let dom_sig_cnt = dom_pct.iter().filter(|p| **p >= 0.1).count();
        v.dom_score_sum = Some(dom_score_sum);
        v.dom_score_max = Some(dom_score_max);
        v.dom_sig_cnt = Some(dom_sig_cnt);
        v.dom_scores = Some(dom);
    }

    // --

    let mut records_by_query = hmms
        .keys()
        .map(|q| (q.as_str(), vec![]))
        .collect::<HashMap<&str, Vec<Record>>>();

    records_by_pair
        .into_iter()
        .for_each(|(k, v)| match records_by_query.get_mut(k.0.as_str()) {
            Some(vec) => vec.push(v),
            None => panic!(),
        });

    records_by_query.values_mut().for_each(|v| {
        v.sort_by(|a, b| match (a.hmmer_p_value, b.hmmer_p_value) {
            (Some(ap), Some(bp)) => ap.partial_cmp(&bp).unwrap(),
            (Some(_), _) => std::cmp::Ordering::Greater,
            (_, Some(_)) => std::cmp::Ordering::Less,
            _ => std::cmp::Ordering::Equal,
        })
    });

    let mut out = BufWriter::new(File::create_new(&args.out_tbl_path)?);

    let widths = HEADER[2]
        .split_whitespace()
        .map(|s| s.len())
        .collect::<Vec<usize>>();

    let mut it = records_by_query.into_values().filter(|v| !v.is_empty());

    let recs = it.next().unwrap();
    write_header(&mut out, &recs)?;

    for recs in it {
        write_records(&mut out, &recs, &widths)?;
    }

    println!(
        "{} took {:.2}s",
        args.out_tbl_path.to_string_lossy(),
        start.elapsed().as_secs_f32()
    );
    Ok(())
}
