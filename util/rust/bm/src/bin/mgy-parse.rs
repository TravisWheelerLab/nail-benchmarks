use std::{
    collections::{HashMap, HashSet},
    fs::{create_dir_all, File, OpenOptions},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{atomic::AtomicUsize, Arc, Mutex},
    time::Instant,
};

use anyhow::Context;
use bioio::tbl::{hmmer::HmmerDomainTable, BlastTable, HitTable, HmmerTable, NailTable};

use clap::{Parser, Subcommand};
use glob::glob;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use regex::Regex;

mod util {
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufRead, BufReader},
        path::Path,
    };

    use anyhow::Context;
    use glob::glob;
    use rayon::ThreadPoolBuilder;

    pub trait Float: PartialOrd {}
    impl Float for f32 {}
    impl Float for f64 {}

    pub fn float_cmp<F: Float>(a: &F, b: &F) -> std::cmp::Ordering {
        a.partial_cmp(b).expect("NaN encountered in float cmp")
    }

    pub fn hit_cmp(a: &bioio::tbl::Hit, b: &bioio::tbl::Hit) -> std::cmp::Ordering {
        float_cmp(&a.e_value, &b.e_value)
    }

    pub fn set_threads(num_threads: usize) -> anyhow::Result<()> {
        ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .context("failed to build rayon global threadpool")
    }

    pub fn score_fmt(f: Option<f32>, w: usize) -> String {
        match f {
            Some(s) => format!("{s:^W$.1}", W = w),
            None => format!("{:^W$}", "-", W = w),
        }
    }

    pub fn p_value_fmt(f: Option<f64>, w: usize) -> String {
        match f {
            Some(s) => format!("{s:^W$.2e}", W = w),
            None => format!("{:^W$}", "-", W = w),
        }
    }

    pub fn int_fmt(f: Option<usize>, w: usize) -> String {
        match f {
            Some(s) => format!("{s:^W$}", W = w),
            None => format!("{:^W$}", "-", W = w),
        }
    }

    pub fn p_value(score: f64, lambda: f64, tau: f64) -> f64 {
        (-lambda * (score - tau)).exp()
    }

    pub struct HmmGumbel {
        pub ga_sc: f32,
        pub ga_p: f64,
        pub tau: f64,
        pub lambda: f64,
    }

    impl HmmGumbel {
        pub fn p_value(&self, score: f64) -> f64 {
            (-self.lambda * (score - self.tau)).exp()
        }
    }

    pub fn target_db_size(target_path: impl AsRef<Path>) -> anyhow::Result<f64> {
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

    pub fn parse_hmms(hmm_path: impl AsRef<Path>) -> anyhow::Result<HashMap<String, HmmGumbel>> {
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
                    .map(|s| s.parse::<f32>().unwrap())
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
                let ga_p = (-lambda * (ga_sc as f64 - tau)).exp();
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

    pub type Cutoffs = HashMap<String, f32>;
    pub fn parse_cutoffs(
        cutoffs_path: impl AsRef<Path>,
        c: usize,
    ) -> anyhow::Result<(Cutoffs, Cutoffs)> {
        let reader = BufReader::new(File::open(cutoffs_path.as_ref())?);

        let mut nail_cutoffs = Cutoffs::new();
        let mut mmseqs_cutoffs = Cutoffs::new();

        for line in reader.lines() {
            let line = line?;

            let (query, rest) = line.split_once(',').unwrap();

            let groups = rest
                .split("),(")
                .map(|g| g.trim_matches(|c| c == '(' || c == ')'))
                .map(|g| {
                    let mut it = g.split(',');
                    let name = it.next().unwrap();
                    let mut nums: Vec<f32> = it.map(|x| x.parse().unwrap()).collect();
                    // note:
                    //   the last number is the
                    //   hit count, not a score
                    nums.pop();
                    (name, nums)
                })
                .collect::<Vec<_>>();

            assert_eq!(groups[0].0, "nail");
            assert_eq!(groups[1].0, "mmseqs");

            let (n, m) = match (groups[0].1.get(c), groups[1].1.get(c)) {
                // only keep cutoffs if:
                //  - there is a cutoff for both, and
                //  - they are both nonzero
                (Some(&n), Some(&m)) if n > 0.0 && m > 0.0 => (n, m),
                _ => continue,
            };

            nail_cutoffs.insert(query.to_string(), n);
            mmseqs_cutoffs.insert(query.to_string(), m);
        }

        Ok((nail_cutoffs, mmseqs_cutoffs))
    }

    pub fn parse_table_indices(dir: impl AsRef<Path>) -> anyhow::Result<Vec<usize>> {
        let dir = dir.as_ref();
        let mut indices = glob(dir.join("*.tbl").to_str().context("invalid *.tbl glob")?)?
            .filter_map(Result::ok)
            .filter_map(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    // NOTE: this assumes the prefixes are of the form
                    //       <tool>.<index>.<blah...>.tbl
                    .and_then(|s| s.split('.').nth(1))
                    .and_then(|i| i.parse::<usize>().ok())
            })
            .collect::<Vec<_>>();

        indices.sort();
        indices.dedup();

        Ok(indices)
    }
}

#[derive(Subcommand)]
enum SubCommands {
    Recall(RecallArgs),
    LearnCutoffs(LearnCutoffsArgs),
    CutoffsSweep(CutoffsSweepArgs),
    Params(ParamsArgs),
}

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    pub command: SubCommands,
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        SubCommands::Recall(args) => recall(args),
        SubCommands::Params(args) => params(args),
        SubCommands::LearnCutoffs(args) => learn_cutoffs(args),
        SubCommands::CutoffsSweep(args) => cutoffs_sweep(args),
    }
}

#[derive(Parser)]
struct CutoffsSweepArgs {
    #[arg(value_name = "nail.tbl")]
    nail_tbl: PathBuf,

    #[arg(value_name = "nail/")]
    nail_dir: PathBuf,

    #[arg(value_name = "mmseqs.tbl")]
    mmseqs_tbl: PathBuf,

    #[arg(value_name = "mmseqs/")]
    mmseqs_dir: PathBuf,

    #[arg(value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(long, value_name = "figures/")]
    figures_dir: Option<PathBuf>,

    #[arg(short = 'e', default_value_t = 1e-3, value_name = "F")]
    reverse_e_cutoff: f64,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,
}

fn cutoffs_sweep(args: CutoffsSweepArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    util::set_threads(args.num_threads)?;

    // ---

    let hmms = util::parse_hmms(&args.query_path)
        .with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    // ---

    let time_pattern = Regex::new(r"Elapsed.*\): (?P<time>.*)").unwrap();

    type DecoyMap = HashMap<String, Vec<bioio::tbl::Hit>>;

    struct Data {
        name: String,
        decoys: DecoyMap,
        saturation: f32,
        time: f32,
    }

    let mut nail_tbl = HitTable::from_path::<_, NailTable>(args.nail_tbl)?;
    nail_tbl.hits.retain(|h| h.e_value <= args.reverse_e_cutoff);
    let nail_map = nail_tbl.to_map();

    let mut nail_data = glob(
        args.nail_dir
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    .map(|path| {
        let mut tbl = HitTable::from_path::<_, NailTable>(&path)
            .unwrap_or_else(|_| panic!("failed to open: {path:?}"));

        let time_path = path.with_extension("time");
        let time: f32 = std::fs::read_to_string(time_path)
            .expect("failed to open .time")
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
            .max_by(util::float_cmp)
            .expect("no times found");

        let name = tbl
            .name
            .strip_suffix(".rev.prf")
            .expect("weird suffix")
            .to_string();

        tbl.hits.retain(|h| h.e_value <= 1.0);

        let mut rev_map = tbl.to_map();
        rev_map.retain(|k, _| !nail_map.contains_key(k));

        let mut decoys = DecoyMap::new();

        rev_map
            .into_values()
            .for_each(|h| decoys.entry(h.query.clone()).or_default().push(h));

        let saturation = decoys.len() as f32 / queries.len() as f32;

        Data {
            name,
            decoys,
            saturation,
            time,
        }
    })
    .collect::<Vec<_>>();

    nail_data.sort_by(|a, b| a.saturation.partial_cmp(&b.saturation).unwrap());
    let w = nail_data.iter().map(|d| d.name.len()).max().unwrap();
    nail_data
        .iter()
        .for_each(|d| println!("{:W$} {:5.1}s {:5.3}", d.name, d.time, d.saturation, W = w));

    // ---

    let mut mmseqs_tbl = HitTable::from_path::<_, BlastTable>(args.mmseqs_tbl)?;
    mmseqs_tbl
        .hits
        .retain(|h| h.e_value <= args.reverse_e_cutoff);
    let mmseqs_map = mmseqs_tbl.to_map();

    let mut mmseqs_data = glob(
        args.mmseqs_dir
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    .map(|path| {
        let tbl = HitTable::from_path::<_, BlastTable>(&path)
            .unwrap_or_else(|_| panic!("failed to open: {path:?}"));

        let time_path = path.with_extension("time");
        let time: f32 = std::fs::read_to_string(time_path)
            .expect("failed to open .time")
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
            .max_by(util::float_cmp)
            .expect("no times found");

        let name = tbl
            .name
            .strip_suffix(".rev.prf")
            .expect("weird suffix")
            .to_string();

        // tbl.hits.retain(|h| h.e_value <= 1.0);

        let mut rev_map = tbl.to_map();
        rev_map.retain(|k, _| !mmseqs_map.contains_key(k));

        let mut decoys = DecoyMap::new();

        rev_map
            .into_values()
            .for_each(|h| decoys.entry(h.query.clone()).or_default().push(h));

        let saturation = decoys.len() as f32 / queries.len() as f32;

        Data {
            name,
            decoys,
            saturation,
            time,
        }
    })
    .collect::<Vec<_>>();

    mmseqs_data.sort_by(|a, b| a.saturation.partial_cmp(&b.saturation).unwrap());
    let w = mmseqs_data.iter().map(|d| d.name.len()).max().unwrap();
    mmseqs_data
        .iter()
        .for_each(|d| println!("{:W$} {:5.1}s {:5.3}", d.name, d.time, d.saturation, W = w));

    println!("took {:.2}s", start.elapsed().as_secs_f32());

    Ok(())
}

#[derive(Parser)]
struct LearnCutoffsArgs {
    #[arg(value_name = "nail/")]
    nail_dir: PathBuf,

    #[arg(value_name = "mmseqs/")]
    mmseqs_dir: PathBuf,

    #[arg(value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(short, long, default_value = "cutoffs.txt", value_name = "cutoffs.txt")]
    out_path: PathBuf,

    #[arg(long, value_name = "figures/")]
    figures_dir: Option<PathBuf>,

    #[arg(short = 'e', default_value_t = 1e-3, value_name = "F")]
    reverse_e_cutoff: f64,

    #[arg(short = 'n', value_name = "N")]
    num_tables: Option<usize>,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,
}

fn learn_cutoffs(args: LearnCutoffsArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    util::set_threads(args.num_threads)?;

    // ---

    let hmms = util::parse_hmms(&args.query_path)
        .with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    let tbl_indices = {
        let indices =
            util::parse_table_indices(&args.nail_dir).context("failed to parse table indices")?;

        if let Some(n) = args.num_tables {
            indices.into_iter().take(n).collect()
        } else {
            indices
        }
    };

    // ---

    type DecoyMap = HashMap<String, Vec<bioio::tbl::Hit>>;

    struct Tables {
        nail: bioio::tbl::HitTable,
        nail_rev: bioio::tbl::HitTable,
        mmseqs: bioio::tbl::HitTable,
        mmseqs_rev: bioio::tbl::HitTable,
    }

    impl Tables {
        fn new(
            nail_dir: impl AsRef<Path>,
            mmseqs_dir: impl AsRef<Path>,
            idx: usize,
        ) -> anyhow::Result<Self> {
            let nail_dir = nail_dir.as_ref();
            let mmseqs_dir = mmseqs_dir.as_ref();

            let path = nail_dir.join(format!("nail.{idx}.prf.tbl"));
            let nail = HitTable::from_path::<_, NailTable>(&path)?;

            let path = nail_dir.join(format!("nail.{idx}.rev.prf.tbl"));
            let nail_rev = HitTable::from_path::<_, NailTable>(&path)?;

            let path = mmseqs_dir.join(format!("mmseqs.{idx}.prf.tbl"));
            let mmseqs = HitTable::from_path::<_, BlastTable>(&path)?;

            let path = mmseqs_dir.join(format!("mmseqs.{idx}.rev.prf.tbl"));
            let mmseqs_rev = HitTable::from_path::<_, BlastTable>(&path)?;

            Ok(Self {
                nail,
                nail_rev,
                mmseqs,
                mmseqs_rev,
            })
        }
    }

    #[derive(Default)]
    struct Data {
        nail_decoys: DecoyMap,
        mmseqs_decoys: DecoyMap,
    }

    impl Data {
        fn update(&mut self, mut tables: Tables, reverse_e_cutoff: f64) {
            // ---
            // filter the real hits by E-value

            tables.nail.hits.retain(|h| h.e_value <= reverse_e_cutoff);
            tables.mmseqs.hits.retain(|h| h.e_value <= reverse_e_cutoff);

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

    let mut data = tbl_indices
        .par_iter()
        .panic_fuse()
        .try_fold(Data::default, |mut data, &idx| -> anyhow::Result<_> {
            let tables = Tables::new(&args.nail_dir, &args.mmseqs_dir, idx)?;
            data.update(tables, args.reverse_e_cutoff);
            Ok(data)
        })
        .try_reduce(Data::default, |mut d1, d2| {
            d1.merge(d2);
            Ok(d1)
        })?;

    // ---

    queries.iter().for_each(|q| {
        data.nail_decoys.entry(q.to_string()).or_default();
        data.mmseqs_decoys.entry(q.to_string()).or_default();
    });

    data.nail_decoys
        .values_mut()
        .for_each(|v| v.sort_by(util::hit_cmp));
    data.mmseqs_decoys
        .values_mut()
        .for_each(|v| v.sort_by(util::hit_cmp));

    // ---

    if let Some(parent) = args.out_path.parent() {
        create_dir_all(parent)?;
    }

    let mut out = BufWriter::new(File::create(&args.out_path)?);

    struct FiguresOut {
        dist: BufWriter<File>,
        ga: BufWriter<File>,
        cnt: BufWriter<File>,
    }
    let mut figures_out = if let Some(figures) = args.figures_dir {
        std::fs::create_dir_all(&figures)?;
        Some(FiguresOut {
            dist: BufWriter::new(File::create(figures.join("dist.txt"))?),
            ga: BufWriter::new(File::create(figures.join("ga.txt"))?),
            cnt: BufWriter::new(File::create(figures.join("count.txt"))?),
        })
    } else {
        None
    };

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
            out,
            "{q},(nail,{},{}),(mmseqs,{},{})",
            nail_scores
                .iter()
                .map(|s| format!("{s:.1}"))
                .collect::<Vec<_>>()
                .join(","),
            nail_hits.len(),
            mmseqs_scores
                .iter()
                .map(|s| format!("{s:.1}"))
                .collect::<Vec<_>>()
                .join(","),
            mmseqs_hits.len(),
        )?;

        // ---

        if let Some(out) = figures_out.as_mut() {
            writeln!(out.dist, "{q},({},{})", nail_hits.len(), mmseqs_hits.len())?;

            writeln!(out.ga, "{:.1},{:.1}", nail_scores[0], hmm.ga_sc)?;

            let diff = hmm.ga_sc - nail_scores[0];
            writeln!(out.cnt, "{},{diff:.1}", nail_hits.len())?;
        }
    }

    // ---

    println!("{:?} took {:?}", args.out_path, start.elapsed());

    Ok(())
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

fn params(args: ParamsArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    let z = util::target_db_size(args.target_path)?;

    let hmms = util::parse_hmms(args.query_path)?;

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
            .max_by(util::float_cmp)
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
            .max_by(util::float_cmp)
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

#[derive(Parser)]
struct RecallArgs {
    #[arg(long, short, value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(long, value_name = "cutoffs.txt")]
    cutoffs_path: PathBuf,

    #[arg(long, value_name = "nail.tbl")]
    nail_path: PathBuf,

    #[arg(long, value_name = "mmseqs.tbl")]
    mmseqs_path: PathBuf,

    #[arg(long, value_name = "hmmer.tbl")]
    hmmer_path: PathBuf,

    #[arg(short, long, value_name = "out/")]
    out_path: PathBuf,

    #[arg(short = 'c', default_value_t = 2usize, value_name = "N")]
    c: usize,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,

    #[arg(short = 'n', value_name = "N")]
    num_tables: Option<usize>,

    #[arg(long)]
    dir: bool,

    #[arg(long)]
    print_times: bool,
}

fn recall(args: RecallArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    util::set_threads(args.num_threads)?;

    // ---

    let (nail_cutoffs, mmseqs_cutoffs) = util::parse_cutoffs(&args.cutoffs_path, args.c)?;

    let hmms = util::parse_hmms(&args.query_path)
        .with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    // ---

    #[derive(Default)]
    struct Times {
        total: AtomicUsize,
        filter: AtomicUsize,
        union: AtomicUsize,
        list: AtomicUsize,
        records: AtomicUsize,
        sort: AtomicUsize,
        read: AtomicUsize,
        write: AtomicUsize,
    }

    impl Times {
        fn print(&self) {
            fn convert(atomic: &AtomicUsize) -> f32 {
                std::time::Duration::from_millis(
                    atomic.load(std::sync::atomic::Ordering::Relaxed) as u64
                )
                .as_secs_f32()
            }

            let total = convert(&self.total);
            let filter = convert(&self.filter);
            let union = convert(&self.union);
            let list = convert(&self.list);
            let records = convert(&self.records);
            let sort = convert(&self.sort);
            let read = convert(&self.read);
            let write = convert(&self.write);
            let misc = total - (filter + union + list + records + sort + read + write);

            println!("filter:  {:5.2}%", filter / total * 100.0);
            println!("union:   {:5.2}%", union / total * 100.0);
            println!("list:    {:5.2}%", list / total * 100.0);
            println!("records: {:5.2}%", records / total * 100.0);
            println!("sort:    {:5.2}%", sort / total * 100.0);
            println!("read:    {:5.2}%", read / total * 100.0);
            println!("write:   {:5.2}%", write / total * 100.0);
            println!("misc:    {:5.2}%", misc / total * 100.0);
        }
    }

    let times: Arc<Times> = Arc::default();

    create_dir_all(&args.out_path)?;

    type Handles = Arc<HashMap<String, Mutex<PathBuf>>>;
    let handles: Handles = Arc::new(
        queries
            .into_iter()
            .map(|query| {
                let path = args.out_path.join(format!("{query}.tbl"));
                let mut file = File::create(&path).unwrap_or_else(|e| {
                    panic!("failed to create output file: {path:?}\n\terror: {e:?}")
                });

                let header = format!(
                    "{:^10}|{:^19}|{:^5}|{:^5}|{:^8}|{:^5}|{:^5}|{:^8}|{:^5}|{:^8}|{:^8}|{:^8}|{:^8}|{}",
                    "query",
                    "target",
                    "n cut",
                    "n sc",
                    "n Eval",
                    "m cut",
                    "m sc",
                    "m Eval",
                    "h sc",
                    "h Eval",
                    "dom max",
                    "dom sum",
                    "dom sig",
                    "dom scores",
                );

                writeln!(file, "{header}").unwrap();
                writeln!(file, "{}", "-".repeat(header.len())).unwrap();

                (query.clone(), Mutex::new(path))
            })
            .collect(),
    );

    // ---

    #[derive(Default)]
    struct Record {
        query: String,
        target: String,
        nail_cutoff: f32,
        nail_score: Option<f32>,
        nail_e_value: Option<f64>,
        mmseqs_cutoff: f32,
        mmseqs_score: Option<f32>,
        mmseqs_e_value: Option<f64>,
        hmmer_score: Option<f32>,
        hmmer_e_value: Option<f64>,
        dom_score_sum: Option<f32>,
        dom_score_max: Option<f32>,
        dom_sig_cnt: Option<usize>,
        dom_scores: Option<Vec<f32>>,
    }

    #[allow(clippy::too_many_arguments)]
    fn process(
        nail_path: impl AsRef<Path>,
        mmseqs_path: impl AsRef<Path>,
        hmmer_path: impl AsRef<Path>,
        nail_cutoffs: &util::Cutoffs,
        mmseqs_cutoffs: &util::Cutoffs,
        handles: Handles,
        times: Arc<Times>,
    ) -> anyhow::Result<()> {
        let start = Instant::now();

        // ---
        let now = Instant::now();

        let nail_map = HitTable::from_path::<_, NailTable>(&nail_path)?.to_map();
        let mmseqs_map = HitTable::from_path::<_, BlastTable>(&mmseqs_path)?.to_map();
        let hmmer_map = HitTable::from_path::<_, HmmerTable>(&hmmer_path)?.to_map();

        times.read.fetch_add(
            now.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        // ---
        let now = Instant::now();

        fn filter_fn(hit: &bioio::tbl::Hit, cutoffs: &util::Cutoffs) -> bool {
            match cutoffs.get(&hit.query) {
                Some(&cutoff) => hit.score > cutoff,
                None => false,
            }
        }

        let mut nail_passed = nail_map
            .iter()
            .filter(|(_, h)| filter_fn(h, nail_cutoffs))
            .map(|(k, _)| k)
            .collect::<HashSet<_>>();

        let mmseqs_passed = mmseqs_map
            .iter()
            .filter(|(_, h)| filter_fn(h, mmseqs_cutoffs))
            .map(|(k, _)| k)
            .collect::<HashSet<_>>();

        let hmmer_passed = hmmer_map
            .iter()
            .filter(|(_, h)| filter_fn(h, nail_cutoffs))
            .map(|(k, _)| k)
            .collect::<HashSet<_>>();

        times.filter.fetch_add(
            now.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        // ---
        let now = Instant::now();

        nail_passed.extend(mmseqs_passed);
        nail_passed.extend(hmmer_passed);
        let passed = nail_passed;

        times.union.fetch_add(
            now.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        // ---
        let now = Instant::now();

        let dom_tbl =
            HmmerDomainTable::from_path(hmmer_path.as_ref().with_extension("domtbl"), |key| {
                passed.contains(&key)
            })?;

        times.read.fetch_add(
            now.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        // ---

        // ---------------------------------------------------------------------
        // if let Some(path) = args.nail_stats_path {
        //     let reader = BufReader::new(File::open(path)?);

        //     let mut it = reader.lines();

        //     while let Ok(batch) = it.by_ref().take(100_000).collect::<Result<Vec<_>, _>>() {
        //         if batch.is_empty() {
        //             break;
        //         }

        //         let filtered = batch
        //             .par_iter()
        //             .filter_map(|line| {
        //                 let mut it = line.split_whitespace();

        //                 let q = it.next()?.to_string();
        //                 let t = it.next()?.to_string();

        //                 // skip 4 fields
        //                 it.nth(3)?;

        //                 let sc = it.next()?.parse::<f32>().ok()?;
        //                 let p = it.next()?.parse::<f64>().ok()?;

        //                 let key = (q, t);
        //                 test.contains(&key).then_some((key, sc, p))
        //             })
        //             .collect::<Vec<_>>();

        //         for (k, sc, p) in filtered {
        //             if let Some(r) = records_by_pair.get_mut(&k) {
        //                 r.nail_cloud_score = Some(sc);
        //                 r.nail_cloud_p_value = Some(p);
        //             }
        //         }
        //     }
        // }
        // ---------------------------------------------------------------------

        let now = Instant::now();

        let mut records_by_query: HashMap<String, Vec<Record>> = HashMap::new();
        for key in passed.into_iter() {
            let query = &key.0;
            let target = &key.1;
            let mut rec = Record {
                query: query.clone(),
                target: target.clone(),
                nail_cutoff: *nail_cutoffs.get(query).expect("no nail cutoff for query"),
                mmseqs_cutoff: *mmseqs_cutoffs
                    .get(query)
                    .expect("no mmseqs cutoff for query"),
                ..Default::default()
            };

            if let Some(h) = nail_map.get(key) {
                rec.nail_score = Some(h.score);
                rec.nail_e_value = Some(h.e_value);
            }

            if let Some(h) = mmseqs_map.get(key) {
                rec.mmseqs_score = Some(h.score);
                rec.mmseqs_e_value = Some(h.e_value);
            }

            if let Some(h) = hmmer_map.get(key) {
                rec.hmmer_score = Some(h.score);
                rec.hmmer_e_value = Some(h.e_value);
            }

            const DOM_SIG_THRESH: f32 = 0.1;
            if let Some(hit) = dom_tbl.hits.get(key) {
                let mut dom_scores = hit.domains.iter().map(|d| d.score).collect::<Vec<_>>();

                let dom_score_sum = dom_scores.iter().sum::<f32>();
                dom_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let dom_score_max = *dom_scores.last().unwrap();

                dom_scores.reverse();

                let dom_pct = dom_scores
                    .iter()
                    .map(|s| s / dom_score_max)
                    .collect::<Vec<_>>();

                let dom_sig_cnt = dom_pct.iter().filter(|p| **p >= DOM_SIG_THRESH).count();

                rec.dom_score_sum = Some(dom_score_sum);
                rec.dom_score_max = Some(dom_score_max);
                rec.dom_sig_cnt = Some(dom_sig_cnt);
                rec.dom_scores = Some(dom_scores);
            }

            records_by_query.entry(query.clone()).or_default().push(rec);
        }

        times.records.fetch_add(
            now.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        // ---
        let now = Instant::now();

        records_by_query.values_mut().for_each(|recs| {
            recs.sort_by(|a, b| {
                if let (Some(sa), Some(sb)) = (a.hmmer_score, b.hmmer_score) {
                    sa.partial_cmp(&sb).unwrap()
                } else if let (Some(_), None) = (a.hmmer_score, b.hmmer_score) {
                    std::cmp::Ordering::Less
                } else if let (None, Some(_)) = (a.hmmer_score, b.hmmer_score) {
                    std::cmp::Ordering::Greater
                } else if let (Some(sa), Some(sb)) = (a.nail_score, b.nail_score) {
                    sa.partial_cmp(&sb).unwrap()
                } else if let (Some(_), None) = (a.nail_score, b.nail_score) {
                    std::cmp::Ordering::Less
                } else if let (None, Some(_)) = (a.nail_score, b.nail_score) {
                    std::cmp::Ordering::Greater
                } else if let (Some(sa), Some(sb)) = (a.mmseqs_score, b.mmseqs_score) {
                    sa.partial_cmp(&sb).unwrap()
                } else if let (Some(_), None) = (a.mmseqs_score, b.mmseqs_score) {
                    std::cmp::Ordering::Less
                } else if let (None, Some(_)) = (a.mmseqs_score, b.mmseqs_score) {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Less
                }
            })
        });

        times.sort.fetch_add(
            now.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        let mut buf: Vec<u8> = vec![];
        let mut remove: Vec<String> = vec![];
        while !records_by_query.is_empty() {
            remove.clear();
            for (query, records) in records_by_query.iter() {
                buf.clear();
                match handles
                    .get(query)
                    .unwrap_or_else(|| panic!("failed to retrive file handle for query: {query}"))
                    .try_lock()
                {
                    Ok(guard) => {
                        let now = Instant::now();

                        let mut file = OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(guard.clone())?;

                        for rec in records {
                            writeln!(
                                buf,
                                "{:10} {:19} {} {} {} {} {} {} {} {} {} {} {} {}",
                                rec.query,
                                rec.target,
                                util::score_fmt(Some(rec.nail_cutoff), 5),
                                util::score_fmt(rec.nail_score, 5),
                                util::p_value_fmt(rec.nail_e_value, 8),
                                util::score_fmt(Some(rec.mmseqs_cutoff), 5),
                                util::score_fmt(rec.mmseqs_score, 5),
                                util::p_value_fmt(rec.mmseqs_e_value, 8),
                                util::score_fmt(rec.hmmer_score, 5),
                                util::p_value_fmt(rec.hmmer_e_value, 8),
                                util::score_fmt(rec.dom_score_max, 5),
                                util::score_fmt(rec.dom_score_sum, 5),
                                util::int_fmt(rec.dom_sig_cnt, 5),
                                match rec.dom_scores {
                                    Some(ref v) => v
                                        .iter()
                                        .map(|f| format!("{f:.1}"))
                                        .collect::<Vec<_>>()
                                        .join(","),
                                    None => "-".to_string(),
                                },
                            )
                            .expect("failed to write record");
                        }

                        file.write_all(&buf)
                            .context("failed to write record buffer to file")?;

                        times.write.fetch_add(
                            now.elapsed().as_millis() as usize,
                            std::sync::atomic::Ordering::Relaxed,
                        );

                        remove.push(query.clone());
                    }
                    Err(_) => continue,
                }
            }

            remove.iter().try_for_each(|q| -> anyhow::Result<()> {
                records_by_query
                    .remove(q)
                    .context("tried to remove a query twice")?;
                Ok(())
            })?;
        }

        times.total.fetch_add(
            start.elapsed().as_millis() as usize,
            std::sync::atomic::Ordering::Relaxed,
        );

        Ok(())
    }

    // ---

    if args.dir {
        let tbl_indices = {
            let indices = util::parse_table_indices(&args.nail_path)
                .context("failed to parse table indices")?;

            if let Some(n) = args.num_tables {
                indices.into_iter().take(n).collect()
            } else {
                indices
            }
        };

        tbl_indices
            .par_iter()
            .panic_fuse()
            .try_for_each(|idx| -> anyhow::Result<()> {
                process(
                    args.nail_path.join(format!("nail.{idx}.prf.tbl")),
                    args.mmseqs_path.join(format!("mmseqs.{idx}.prf.tbl")),
                    args.hmmer_path.join(format!("hmmer.{idx}.prf.tbl")),
                    &nail_cutoffs,
                    &mmseqs_cutoffs,
                    handles.clone(),
                    times.clone(),
                )?;
                Ok(())
            })?;
    } else {
        process(
            args.nail_path,
            args.mmseqs_path,
            args.hmmer_path,
            &nail_cutoffs,
            &mmseqs_cutoffs,
            handles,
            times.clone(),
        )?;
    }

    println!(
        "{} took {:.2}s",
        args.out_path.to_string_lossy(),
        start.elapsed().as_secs_f32()
    );

    if args.print_times {
        times.print();
    }

    Ok(())
}
