use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    fs::{create_dir_all, File, OpenOptions},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{atomic::AtomicUsize, Arc, Mutex},
    time::Instant,
};

use anyhow::{bail, Context};
use bioio::tbl::{hmmer::HmmerDomainTable, BlastTable, HitTable, HmmerTable, NailTable};

use clap::{Parser, Subcommand};
use glob::glob;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use regex::Regex;
use thread_local::ThreadLocal;

mod util {
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufRead, BufReader},
        path::{Path, PathBuf},
        sync::{Arc, Mutex},
    };

    use anyhow::Context;
    use bioio::tbl::{BlastTable, HitTable, HmmerTable, NailTable};
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

    pub fn parse_table_indices(
        dir: impl AsRef<Path>,
        n: Option<usize>,
    ) -> anyhow::Result<Vec<usize>> {
        let dir = dir.as_ref();
        let mut indices = glob(dir.join("*.tbl").to_str().context("invalid *.tbl glob")?)?
            .filter_map(Result::ok)
            // the runs table shares the extension but is not a hit table
            .filter(|p| p.file_name().is_some_and(|n| n != run::table::FILE_NAME))
            .filter_map(|p| {
                // the shard index is the last dot-separated token of the stem.
                // It cannot be taken from the front: run names embed floating
                // point parameters, so `nail-s12.0-prog.prf.7` has dots well
                // before the index.
                p.file_stem()
                    .and_then(|n| n.to_str())
                    .and_then(|s| s.rsplit('.').next())
                    .and_then(|i| i.parse::<usize>().ok())
            })
            .collect::<Vec<_>>();

        indices.sort();
        indices.dedup();

        if indices.is_empty() {
            anyhow::bail!(
                "found no shard-indexed hit tables in {}; \
                 expected names like <run>.<shard>.tbl",
                dir.display()
            );
        }

        if let Some(n) = n {
            Ok(indices.into_iter().take(n).collect())
        } else {
            Ok(indices)
        }
    }

    pub struct Tables {
        pub nail: bioio::tbl::HitTable,
        pub nail_rev: bioio::tbl::HitTable,
        pub mmseqs: bioio::tbl::HitTable,
        pub mmseqs_rev: bioio::tbl::HitTable,
        pub hmmer: bioio::tbl::HitTable,
        pub hmmer_rev: bioio::tbl::HitTable,
    }

    impl Tables {
        pub fn new(
            nail_dir: impl AsRef<Path>,
            mmseqs_dir: impl AsRef<Path>,
            hmmer_dir: impl AsRef<Path>,
            query: &str,
        ) -> anyhow::Result<Self> {
            let nail_dir = nail_dir.as_ref();
            let mmseqs_dir = mmseqs_dir.as_ref();
            let hmmer_dir = hmmer_dir.as_ref();

            let path = nail_dir.join(format!("{query}.tbl"));
            let nail = HitTable::from_path::<_, NailTable>(&path)?;

            let path = nail_dir.join(format!("{query}.rev.tbl"));
            let nail_rev = HitTable::from_path::<_, NailTable>(&path)?;

            let path = mmseqs_dir.join(format!("{query}.tbl"));
            let mmseqs = HitTable::from_path::<_, BlastTable>(&path)?;

            let path = mmseqs_dir.join(format!("{query}.rev.tbl"));
            let mmseqs_rev = HitTable::from_path::<_, BlastTable>(&path)?;

            let path = hmmer_dir.join(format!("{query}.tbl"));
            let hmmer = HitTable::from_path::<_, HmmerTable>(&path)?;

            let path = hmmer_dir.join(format!("{query}.rev.tbl"));
            let hmmer_rev = HitTable::from_path::<_, HmmerTable>(&path)?;

            Ok(Self {
                nail,
                nail_rev,
                mmseqs,
                mmseqs_rev,
                hmmer,
                hmmer_rev,
            })
        }
    }

    #[derive(Clone)]
    pub struct PathHandles {
        handles: Arc<HashMap<String, Mutex<PathBuf>>>,
    }

    impl PathHandles {
        pub fn new<I, S, P, F>(keys: I, dir: P, ext: &str, mut init: F) -> Self
        where
            I: IntoIterator<Item = S>,
            S: AsRef<str>,
            P: AsRef<Path>,
            F: FnMut(&PathBuf),
        {
            let dir = dir.as_ref();
            Self {
                handles: Arc::new(
                    keys.into_iter()
                        .map(|name: S| {
                            let name = name.as_ref();
                            let path = dir.join(name).with_extension(ext);
                            init(&path);

                            (name.to_string(), Mutex::new(path))
                        })
                        .collect(),
                ),
            }
        }

        pub fn get(&self, key: &str) -> Option<&Mutex<PathBuf>> {
            self.handles.get(key)
        }
    }
}

/// Analysis subcommands for this benchmark. Kept here rather than in the `run`
/// library because what counts as a result is benchmark-specific.
#[derive(Subcommand)]
pub enum Cmd {
    Recall(RecallArgs),
    LearnCutoffs(LearnCutoffsArgs),
    CutoffsSweep(CutoffsSweepArgs),
    Params(ParamsArgs),
    CheckRev(CheckRevArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Recall(args) => recall(args),
        Cmd::Params(args) => params(args),
        Cmd::LearnCutoffs(args) => learn_cutoffs(args),
        Cmd::CutoffsSweep(args) => cutoffs_sweep(args),
        Cmd::CheckRev(args) => check_rev(args),
    }
}

#[derive(Parser)]
pub struct CutoffsSweepArgs {
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

    let hmms = bioio::hmm::parse_stats(&args.query_path)
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
pub struct LearnCutoffsArgs {
    #[arg(value_name = "nail/")]
    nail_dir: PathBuf,

    #[arg(value_name = "mmseqs/")]
    mmseqs_dir: PathBuf,

    #[arg(value_name = "hmmer/")]
    hmmer_dir: PathBuf,

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

    let hmms = bioio::hmm::parse_stats(&args.query_path)
        .with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    if let Some(parent) = args.out_path.parent() {
        create_dir_all(parent)?;
    }

    let out = Arc::new(Mutex::new(BufWriter::new(File::create(&args.out_path)?)));

    queries.par_iter().try_for_each(|q| -> anyhow::Result<()> {
        let mut tables =
            match util::Tables::new(&args.nail_dir, &args.mmseqs_dir, &args.hmmer_dir, q) {
                Ok(t) => t,
                Err(_) => return Ok(()),
            };

        // ---
        // filter the real hits by E-value
        tables
            .nail
            .hits
            .retain(|h| h.e_value <= args.reverse_e_cutoff);

        tables
            .mmseqs
            .hits
            .retain(|h| h.e_value <= args.reverse_e_cutoff);

        tables
            .hmmer
            .hits
            .retain(|h| h.e_value <= args.reverse_e_cutoff);

        // ---
        // convert to (query, target)-keyed maps for easier comparison

        let nail_map = tables.nail.to_map();
        let mut nail_rev_map = tables.nail_rev.to_map();

        let mmseqs_map = tables.mmseqs.to_map();
        let mut mmseqs_rev_map = tables.mmseqs_rev.to_map();

        let hmmer_map = tables.hmmer.to_map();
        let mut hmmer_rev_map = tables.hmmer_rev.to_map();

        // ---
        // retain only reverse hits for pairs that don't
        // remain in the real hits after filtering

        nail_rev_map.retain(|k, _| !nail_map.contains_key(k));
        mmseqs_rev_map.retain(|k, _| !mmseqs_map.contains_key(k));
        hmmer_rev_map.retain(|k, _| !hmmer_map.contains_key(k));

        let mut n = nail_rev_map.into_values().collect::<Vec<_>>();
        let mut m = mmseqs_rev_map.into_values().collect::<Vec<_>>();
        let mut h = hmmer_rev_map.into_values().collect::<Vec<_>>();

        n.sort_by(util::hit_cmp);
        m.sort_by(util::hit_cmp);
        h.sort_by(util::hit_cmp);

        const N_CUTOFF: usize = 5;

        let nn = n
            .iter()
            .map(|t| t.score)
            .chain(std::iter::repeat(0.0))
            .take(N_CUTOFF)
            .collect::<Vec<_>>();

        let mm = m
            .iter()
            .map(|t| t.score)
            .chain(std::iter::repeat(0.0))
            .take(N_CUTOFF)
            .collect::<Vec<_>>();

        let hh = h
            .iter()
            .map(|t| t.score)
            .chain(std::iter::repeat(0.0))
            .take(N_CUTOFF)
            .collect::<Vec<_>>();

        match out.lock() {
            Ok(mut guard) => {
                writeln!(
                    guard,
                    "{q},(nail,{},{}),(mmseqs,{},{})(hmmer,{},{})",
                    nn.iter()
                        .map(|s| format!("{s:.1}"))
                        .collect::<Vec<_>>()
                        .join(","),
                    n.len(),
                    mm.iter()
                        .map(|s| format!("{s:.1}"))
                        .collect::<Vec<_>>()
                        .join(","),
                    m.len(),
                    hh.iter()
                        .map(|s| format!("{s:.1}"))
                        .collect::<Vec<_>>()
                        .join(","),
                    h.len(),
                )?;
            }
            Err(_) => panic!("poisoned"),
        }

        Ok(())
    })?;

    println!("{:?} took {:?}", args.out_path, start.elapsed());

    Ok(())
}

#[derive(Parser)]
pub struct ParamsArgs {
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

    #[arg(long, value_name = "hmmer.domtbl")]
    hmmer_dom_path: PathBuf,

    #[arg(short, long, value_name = "out/")]
    out_path: PathBuf,

    #[arg(short = 'c', default_value_t = 2usize, value_name = "N")]
    c: usize,

    #[arg(long)]
    print_times: bool,
}

fn params(args: ParamsArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    // util::set_threads(args.num_threads)?;

    // ---

    let (nail_cutoffs, mmseqs_cutoffs) =
        util::parse_cutoffs(&args.cutoffs_path, args.c).context("cutoffs")?;

    let hmms = bioio::hmm::parse_stats(&args.query_path)
        .with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    // ---

    let hmmer_tbl = HitTable::from_path::<_, HmmerTable>(args.hmmer_path)?;
    let hmmer_dom = HmmerDomainTable::from_path(args.hmmer_dom_path, |_| true)?;

    let mut passed_hmmer = HashSet::new();
    let mut hmmer_cnt = 0.0;
    let mut fuck = 0;
    for hit in hmmer_tbl.hits {
        let cutoff = match nail_cutoffs.get(&hit.query) {
            Some(c) => c,
            None => continue,
        };

        let domain_hit = match hmmer_dom.hits.get(&(hit.query.clone(), hit.target.clone())) {
            Some(h) => h,
            None => {
                fuck += 1;
                continue;
            }
        };

        if domain_hit.domains.iter().any(|h| h.score >= *cutoff) {
            passed_hmmer.insert((hit.query, hit.target));
            hmmer_cnt += 1.0;
        }
    }

    println!("{fuck}");

    // ---

    let time_pattern = Regex::new(r"Elapsed.*\): (?P<time>.*)$").unwrap();

    let mut out = BufWriter::new(File::create(&args.out_path)?);

    for path in glob(
        args.nail_path
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    {
        let mut cnt = 0.0;

        let tbl = HitTable::from_path::<_, NailTable>(&path)?;
        for hit in tbl.hits {
            let cutoff = match nail_cutoffs.get(&hit.query) {
                Some(c) => c,
                None => continue,
            };

            if passed_hmmer.contains(&(hit.query, hit.target)) && hit.score >= *cutoff {
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

        let name = format!("{prefix}.{search_type}");

        writeln!(out, "{name},({time:.4},{frac:.4})")?;
    }

    for path in glob(
        args.mmseqs_path
            .join("*.tbl")
            .to_str()
            .context("invalid *.tbl glob")?,
    )?
    .filter_map(Result::ok)
    {
        let mut cnt = 0.0;

        let tbl = HitTable::from_path::<_, BlastTable>(&path)?;
        for hit in tbl.hits {
            let cutoff = match mmseqs_cutoffs.get(&hit.query) {
                Some(c) => c,
                None => continue,
            };
            if passed_hmmer.contains(&(hit.query, hit.target)) && hit.score >= *cutoff {
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
                    .map(|t| t.parse::<f32>().expect("failed to parse time"))
                    .rev()
                    .enumerate()
                    .map(|(i, t)| t * 60.0_f32.powf(i as f32))
                    .sum()
            })
            .max_by(util::float_cmp)
            .expect("no times found");

        let name = format!("{prefix}.{search_type}");
        writeln!(out, "{name},({time:.4},{frac:.4})")?;
    }

    println!(
        "{} took {:.2}s",
        args.out_path.to_string_lossy(),
        start.elapsed().as_secs_f32()
    );

    Ok(())
}

#[derive(Parser)]
pub struct RecallArgs {
    #[arg(long, short, value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(long, value_name = "cutoffs.txt")]
    cutoffs_path: PathBuf,

    /// Which benchmark directory under benchmarks/mgnify/ to read.
    #[arg(long, default_value = crate::build::DEFAULT_NAME)]
    name: String,

    /// Directory holding the per-shard hit tables.
    #[arg(long, value_name = "results/")]
    results: Option<PathBuf>,

    /// The runs table indexing those results.
    #[arg(long, value_name = "runs.tbl")]
    runs: Option<PathBuf>,

    /// Run name to read for each tool. Only needed when the config sweeps a
    /// tool over several parameter settings; otherwise the single run for that
    /// tool is used.
    #[arg(long, value_name = "RUN")]
    nail: Option<String>,

    #[arg(long, value_name = "RUN")]
    mmseqs: Option<String>,

    #[arg(long, value_name = "RUN")]
    hmmer: Option<String>,

    #[arg(short, long, value_name = "out/")]
    out_path: Option<PathBuf>,

    #[arg(short = 'c', default_value_t = 2usize, value_name = "N")]
    c: usize,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,

    /// Only read this many shards.
    #[arg(short = 'n', value_name = "N")]
    num_tables: Option<usize>,

    #[arg(long)]
    print_times: bool,
}

fn recall(args: RecallArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    // a benchmark owns its results, runs table, and analysis output
    let bench = crate::build::dir().join(&args.name);
    let results_dir = args.results.unwrap_or_else(|| bench.join("results"));
    let runs_path = args.runs.unwrap_or_else(|| bench.join(run::table::FILE_NAME));
    let out_path = args.out_path.unwrap_or_else(|| bench.join("analysis"));

    util::set_threads(args.num_threads)?;

    // ---

    let (nail_cutoffs, mmseqs_cutoffs) = util::parse_cutoffs(&args.cutoffs_path, args.c)?;

    let hmms = bioio::hmm::parse_stats(&args.query_path)
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

    create_dir_all(&out_path)?;

    let header = format!(
        "{:^10}|{:^19}|{:^8}|{:^5}|{:^5}|{:^8}|{:^5}|{:^5}|{:^8}|{:^5}|{:^5}|{:^8}|{:^8}|{:^8}|{:^8}|{}",
        "query",
        "target",
        "file",
        "n cut",
        "n sc",
        "n Eval",
        "m cut",
        "m sc",
        "m Eval",
        "h cut",
        "h sc",
        "h Eval",
        "dom max",
        "dom sum",
        "dom sig",
        "dom scores",
    );

    let handles = util::PathHandles::new(queries, &out_path, "tbl", |p| {
        let mut f = File::create(p)
            .unwrap_or_else(|e| panic!("failed to create output file: {p:?}\n\terror: {e:?}"));
        writeln!(f, "{header}").unwrap();
        writeln!(f, "{}", "-".repeat(header.len())).unwrap();

        let mut f = File::create(p.with_extension("md.tbl"))
            .unwrap_or_else(|e| panic!("failed to create output file: {p:?}\n\terror: {e:?}"));
        writeln!(f, "{header}").unwrap();
        writeln!(f, "{}", "-".repeat(header.len())).unwrap();
    });

    // ---

    struct Stats {
        map: HashMap<String, AtomicUsize>,
    }

    impl Stats {
        fn new() -> Self {
            Self {
                map: [
                    ("nail".to_string(), AtomicUsize::new(0)),
                    ("nail_hmmer".to_string(), AtomicUsize::new(0)),
                    ("nail_hmmer_single".to_string(), AtomicUsize::new(0)),
                    ("mmseqs".to_string(), AtomicUsize::new(0)),
                    ("mmseqs_hmmer".to_string(), AtomicUsize::new(0)),
                    ("mmseqs_hmmer_single".to_string(), AtomicUsize::new(0)),
                    ("hmmer".to_string(), AtomicUsize::new(0)),
                    ("hmmer_single".to_string(), AtomicUsize::new(0)),
                ]
                .into_iter()
                .collect(),
            }
        }

        fn add(&self, key: &str, val: usize) {
            self.map
                .get(key)
                .unwrap_or_else(|| panic!("no stats key: {key}"))
                .fetch_add(val, std::sync::atomic::Ordering::Relaxed);
        }

        fn get(&self, key: &str) -> f32 {
            let atomic = self.map.get(key).unwrap();
            atomic.load(std::sync::atomic::Ordering::SeqCst) as f32
        }

        fn print(&self) {
            let nail = self.get("nail");
            let nail_hmmer = self.get("nail_hmmer");
            let nail_hmmer_single = self.get("nail_hmmer_single");
            let mmseqs = self.get("mmseqs");
            let mmseqs_hmmer = self.get("mmseqs_hmmer");
            let mmseqs_hmmer_single = self.get("mmseqs_hmmer_single");
            let hmmer = self.get("hmmer");
            let hmmer_single = self.get("hmmer_single");

            println!("total:");
            println!("hmmer:  {hmmer} 1.0");
            println!("nail:   {nail} {:.3}", nail_hmmer / hmmer);
            println!("mmseqs: {mmseqs} {:.3}", mmseqs_hmmer / hmmer);

            println!();
            println!("single domain:");
            println!("hmmer:  {hmmer_single} 1.0");
            println!(
                "nail:   {nail_hmmer_single} {:.3}",
                nail_hmmer_single / hmmer_single
            );
            println!(
                "mmseqs: {mmseqs_hmmer_single} {:.3}",
                mmseqs_hmmer_single / hmmer_single
            );
        }
    }

    let stats = Arc::new(Stats::new());

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
        hmmer_cutoff: f32,
        hmmer_score: Option<f32>,
        hmmer_e_value: Option<f64>,
        dom_score_sum: Option<f32>,
        dom_score_max: Option<f32>,
        dom_sig_cnt: Option<usize>,
        dom_scores: Option<Vec<f32>>,
        file: String,
    }

    impl Record {
        pub fn write(&self, buf: &mut impl Write) {
            writeln!(
                buf,
                "{:10} {:19} {:8} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                self.query,
                self.target,
                self.file,
                util::score_fmt(Some(self.nail_cutoff), 5),
                util::score_fmt(self.nail_score, 5),
                util::p_value_fmt(self.nail_e_value, 8),
                util::score_fmt(Some(self.mmseqs_cutoff), 5),
                util::score_fmt(self.mmseqs_score, 5),
                util::p_value_fmt(self.mmseqs_e_value, 8),
                util::score_fmt(Some(self.hmmer_cutoff), 5),
                util::score_fmt(self.hmmer_score, 5),
                util::p_value_fmt(self.hmmer_e_value, 8),
                util::score_fmt(self.dom_score_max, 5),
                util::score_fmt(self.dom_score_sum, 5),
                util::int_fmt(self.dom_sig_cnt, 5),
                match self.dom_scores {
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
    }

    #[allow(clippy::too_many_arguments)]
    fn process(
        // the target file name, as it appears in the runs table; deriving it
        // from the hit table's path is what produced "0-prog" out of
        // nail-s12.0-prog.prf.1.tbl
        shard: &str,
        nail_path: impl AsRef<Path>,
        mmseqs_path: impl AsRef<Path>,
        hmmer_path: impl AsRef<Path>,
        nail_cutoffs: &util::Cutoffs,
        mmseqs_cutoffs: &util::Cutoffs,
        handles: util::PathHandles,
        stats: Arc<Stats>,
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

        let now = Instant::now();

        let mut nail = 0;
        let mut nail_hmmer = 0;
        let mut nail_hmmer_single = 0;
        let mut mmseqs = 0;
        let mut mmseqs_hmmer = 0;
        let mut mmseqs_hmmer_single = 0;
        let mut hmmer = 0;
        let mut hmmer_single = 0;

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
                // shard is already the target file name from the runs table
                file: shard.to_string(),
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

            let n = rec.nail_score.unwrap_or(0.0);
            let m = rec.mmseqs_score.unwrap_or(0.0);
            let h = rec.hmmer_score.unwrap_or(0.0);
            let s = rec.dom_sig_cnt.unwrap_or(0);

            if n > rec.nail_cutoff {
                nail += 1;
            }

            if m > rec.mmseqs_cutoff {
                mmseqs += 1;
            }

            if h > rec.nail_cutoff {
                hmmer += 1;
                if s == 1 {
                    hmmer_single += 1;
                }
            }

            if n > rec.nail_cutoff && h > rec.nail_cutoff {
                nail_hmmer += 1;
                if s == 1 {
                    nail_hmmer_single += 1;
                }
            }

            if m > rec.mmseqs_cutoff && h > rec.nail_cutoff {
                mmseqs_hmmer += 1;
                if s == 1 {
                    mmseqs_hmmer_single += 1;
                }
            }

            records_by_query.entry(query.clone()).or_default().push(rec);
        }

        stats.add("nail", nail);
        stats.add("mmseqs", mmseqs);
        stats.add("hmmer", hmmer);

        stats.add("nail_hmmer", nail_hmmer);
        stats.add("mmseqs_hmmer", mmseqs_hmmer);

        stats.add("nail_hmmer_single", nail_hmmer_single);
        stats.add("mmseqs_hmmer_single", mmseqs_hmmer_single);

        stats.add("hmmer_single", hmmer_single);

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
        let mut md_buf: Vec<u8> = vec![];
        let mut remove: Vec<String> = vec![];
        while !records_by_query.is_empty() {
            remove.clear();
            for (query, records) in records_by_query.iter() {
                buf.clear();
                md_buf.clear();
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

                        let mut md_file = OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(guard.with_extension("md.tbl"))?;

                        for rec in records {
                            rec.write(&mut buf);

                            if rec.dom_sig_cnt >= Some(2) && rec.nail_score == None {
                                rec.write(&mut md_buf);
                            }
                        }

                        file.write_all(&buf)
                            .context("failed to write record buffer to file")?;

                        md_file
                            .write_all(&md_buf)
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

    // the runs table is the index of a results directory: it records which run
    // names exist and which shards each covered, so nothing has to be inferred
    // from filenames. It sits beside benchmark/ rather than inside results/,
    // so the two paths are given separately.
    let runs = run::Runs::load(&runs_path, &results_dir)?;

    let pick = |explicit: Option<String>, tool: &str| -> anyhow::Result<String> {
        match explicit {
            Some(name) => {
                if runs.targets(&name).is_empty() {
                    anyhow::bail!(
                        "no run named {name:?} in {}; available: {}",
                        runs_path.display(),
                        runs.names().join(", ")
                    );
                }
                Ok(name)
            }
            None => Ok(runs.only_for_tool(tool)?.to_string()),
        }
    };

    let nail_run = pick(args.nail, "nail")?;
    let mmseqs_run = pick(args.mmseqs, "mmseqs")?;
    let hmmer_run = pick(args.hmmer, "hmmer")?;

    println!("nail:   {nail_run}");
    println!("mmseqs: {mmseqs_run}");
    println!("hmmer:  {hmmer_run}");

    // only shards all three tools completed are comparable
    let mut shards = runs.shared_targets(&[&nail_run, &mmseqs_run, &hmmer_run]);
    if shards.is_empty() {
        anyhow::bail!("no shard was completed by all three of nail, mmseqs and hmmer");
    }
    if let Some(n) = args.num_tables {
        shards.truncate(n);
    }
    println!("{} shards\n", shards.len());

    shards
        .par_iter()
        .panic_fuse()
        .try_for_each(|shard| -> anyhow::Result<()> {
            process(
                shard,
                runs.table_path(&nail_run, shard),
                runs.table_path(&mmseqs_run, shard),
                runs.table_path(&hmmer_run, shard),
                &nail_cutoffs,
                &mmseqs_cutoffs,
                handles.clone(),
                stats.clone(),
                times.clone(),
            )?;
            Ok(())
        })?;

    println!(
        "{} took {:.2}s",
        out_path.to_string_lossy(),
        start.elapsed().as_secs_f32()
    );

    println!();

    if args.print_times {
        times.print();
    }

    stats.print();

    Ok(())
}

#[derive(Parser)]
pub struct CheckRevArgs {
    #[arg(long, value_name = "nail/")]
    nail_dir: PathBuf,

    #[arg(long, value_name = "mmseqs/")]
    mmseqs_dir: PathBuf,

    #[arg(long, value_name = "query.hmm")]
    query_path: PathBuf,

    #[arg(long, value_name = "targets/")]
    target_dir: PathBuf,

    #[arg(short, long, default_value = "fwd-rev/", value_name = "fwd-rev/")]
    out_path: PathBuf,

    #[arg(short = 'n', value_name = "N")]
    num_tables: Option<usize>,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,
}

fn check_rev(args: CheckRevArgs) -> anyhow::Result<()> {
    let start = Instant::now();

    util::set_threads(args.num_threads)?;

    // ---

    let hmms = bioio::hmm::parse_stats(&args.query_path)
        .with_context(|| format!("failed to open: {:?}", args.query_path))?;

    let mut queries = hmms.keys().collect::<Vec<_>>();
    queries.sort();

    create_dir_all(&args.out_path)?;

    let handles = util::PathHandles::new(queries, &args.out_path, "fa", |p| {
        File::create(p)
            .unwrap_or_else(|e| panic!("failed to create output file: {p:?}\n\terror: {e:?}"));
    });

    // ---

    let tl_buf: ThreadLocal<RefCell<Vec<u8>>> = ThreadLocal::new();

    let tbl_indices = util::parse_table_indices(&args.nail_dir, args.num_tables)?;
    tbl_indices
        .par_iter()
        .panic_fuse()
        .try_for_each(|&idx| -> anyhow::Result<()> {
            let target_path = args.target_dir.join(format!("{idx}.fa"));
            let target = bioio::fasta::Fasta::from_path(target_path)?;

            // ---

            let nail_tbl = bioio::tbl::HitTable::from_path::<_, NailTable>(
                args.nail_dir.join(format!("nail.{idx}.rev.prf.tbl")),
            )?;

            let mmseqs_tbl = bioio::tbl::HitTable::from_path::<_, BlastTable>(
                args.mmseqs_dir.join(format!("mmseqs.{idx}.rev.prf.tbl")),
            )?;

            // ---

            fn map_fn(tbl: bioio::tbl::HitTable, map: &mut HashMap<String, Vec<String>>) {
                tbl.to_query_map().into_iter().for_each(|(q, hits)| {
                    let v = map.entry(q).or_default();
                    hits.iter().for_each(|h| v.push(h.target.clone()));
                });
            }

            let mut map: HashMap<String, Vec<String>> = HashMap::new();

            map_fn(nail_tbl, &mut map);
            map_fn(mmseqs_tbl, &mut map);

            map.values_mut().for_each(|targets| {
                targets.sort();
                targets.dedup();
            });

            // ---

            let mut buf = tl_buf.get_or(|| RefCell::new(vec![])).borrow_mut();

            map.into_iter()
                .try_for_each(|(query, targets)| -> anyhow::Result<()> {
                    targets.iter().try_for_each(|t| -> anyhow::Result<()> {
                        let seq = target.records.get(t).unwrap();
                        writeln!(buf, "{seq}")?;
                        Ok(())
                    })?;

                    match handles.get(&query).context("")?.lock() {
                        Ok(path) => {
                            let mut file = OpenOptions::new()
                                .create(true)
                                .append(true)
                                .open(path.clone())?;

                            file.write_all(&buf)?;
                            buf.clear();
                        }
                        Err(_) => bail!("mutex poisoned"),
                    }

                    Ok(())
                })?;

            Ok(())
        })?;

    println!(
        "{} took {:.2}s",
        args.out_path.to_string_lossy(),
        start.elapsed().as_secs_f32()
    );

    Ok(())
}
