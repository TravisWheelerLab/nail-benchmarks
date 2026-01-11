use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{stdout, BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

use anyhow::Context;
use clap::Parser;
use rayon::prelude::*;

use crate::mmseqs_db::{Db, SplitDb};

#[derive(Parser)]
struct Args {
    #[arg(value_name = "dir/")]
    dir_path: PathBuf,

    #[arg(value_name = "ga_hits.tbl")]
    ga_hits: PathBuf,

    #[arg(short, long, value_name = "path")]
    out_path: Option<PathBuf>,

    #[arg(short, value_name = "n", default_value_t = 8)]
    threads: usize,
}

mod mmseqs_db {
    use std::{
        collections::HashMap,
        fs::{self, File},
        io::{BufRead, BufReader, Read, Seek},
        path::{Path, PathBuf},
    };

    use anyhow::Context;
    use regex::Regex;

    pub struct Index {
        pub offsets: Vec<u64>,
        pub lengths: Vec<u64>,
    }

    impl Index {
        pub fn from_path<P>(path: P) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
        {
            let reader =
                BufReader::new(File::open(&path).with_context(|| format!("{:?}", path.as_ref()))?);

            let mut offsets = vec![];
            let mut lengths = vec![];
            for line in reader.lines() {
                let line = line?;

                let tokens = line
                    .split('\t')
                    .map(|s| s.parse::<u64>().expect("failed to parse index line"))
                    .collect::<Vec<_>>();

                offsets.push(tokens[1]);
                // -1 because the index lengths include a null byte
                lengths.push(tokens[2] - 1);
            }

            offsets.shrink_to_fit();
            lengths.shrink_to_fit();

            Ok(Self { offsets, lengths })
        }
    }

    pub struct Db {
        pub index: Index,
        pub path: PathBuf,
    }

    impl Db {
        pub fn from_path<P>(path: P) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
        {
            let path = path.as_ref();
            let index_path = path.with_extension("index");

            Ok(Self {
                index: Index::from_path(index_path)?,
                path: path.to_path_buf(),
            })
        }

        pub fn len(&self) -> usize {
            self.index.offsets.len()
        }

        pub fn get(&self, idx: usize) -> anyhow::Result<String> {
            let offset = self.index.offsets[idx];
            let length = self.index.lengths[idx];

            let mut file = File::open(&self.path).with_context(|| format!("{:?}", &self.path))?;
            file.seek(std::io::SeekFrom::Start(offset))?;
            let mut buf = vec![0u8; length as usize];
            file.read_exact(&mut buf)?;

            let mut value = String::from_utf8(buf)?;
            value.truncate(value.trim_end().len());
            Ok(value)
        }

        pub fn into_map(self) -> HashMap<String, usize> {
            let mut map = HashMap::new();
            for idx in 0..self.len() {
                let key = self.get(idx).unwrap();
                map.insert(key, idx);
            }

            map
        }
    }

    pub struct SplitDb {
        pub index: Index,
        pub paths: Vec<PathBuf>,
        pub sizes: Vec<u64>,
    }

    impl SplitDb {
        pub fn from_path<P>(path: P) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
        {
            let path = path.as_ref();
            let index_path = path.with_extension("index");

            let dir = path.parent().unwrap();
            let name = path.file_stem().unwrap().to_str().unwrap();

            let re = Regex::new(&format!(r"^{}\.\d+$", regex::escape(name)))?;
            let mut paths = fs::read_dir(dir)?
                .filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    p.file_name()
                        .and_then(|s| s.to_str())
                        .map(|s| re.is_match(s))
                        .unwrap_or(false)
                })
                .collect::<Vec<_>>();

            paths.sort_by_key(|p| {
                p.extension()
                    .and_then(|e| e.to_str())
                    .and_then(|s| s.parse::<usize>().ok())
                    .expect("bad DB split suffix")
            });

            let sizes = paths
                .iter()
                .map(|p| fs::metadata(p).unwrap().len())
                .collect::<Vec<_>>();

            let prefix = sizes
                .into_iter()
                .scan(0u64, |acc, x| {
                    let cur = *acc;
                    *acc += x;
                    Some(cur)
                })
                .collect::<Vec<_>>();

            Ok(Self {
                index: Index::from_path(index_path)?,
                paths,
                sizes: prefix,
            })
        }

        pub fn len(&self) -> usize {
            self.index.offsets.len()
        }

        fn offset_idx(&self, offset: u64) -> usize {
            match self.sizes.binary_search(&offset) {
                Ok(idx) => idx,
                Err(idx) => idx,
            }
        }

        pub fn get(&self, idx: usize) -> anyhow::Result<String> {
            let offset = self.index.offsets[idx];
            let length = self.index.lengths[idx];

            let file_idx = match self.sizes.binary_search(&offset) {
                // this means the binary search found the exact offset
                Ok(idx) => idx,
                // this means the binary search landed between offsets
                Err(idx) => idx.saturating_sub(1),
            };

            let offset = offset - self.sizes[file_idx];
            let path = &self.paths[file_idx];

            let mut file = File::open(path).with_context(|| format!("{path:?}"))?;
            file.seek(std::io::SeekFrom::Start(offset))?;
            let mut buf = vec![0u8; length as usize];
            file.read_exact(&mut buf)?;

            Ok(String::from_utf8(buf)?)
        }
    }
}

pub fn print_histogram(data: &[(usize, usize)]) {
    if data.is_empty() {
        return;
    }

    let max = data.iter().map(|&(_, v)| v).max().unwrap();
    let maxw = max.to_string().len();
    let width = 40usize;
    let scale = max.div_ceil(width);
    let scale = scale.max(1);

    for &(label, value) in data {
        let blocks = value / scale;
        println!(
            "{:>4}b | {value:>W$} | {}",
            label,
            "█".repeat(blocks),
            W = maxw
        );
    }
}

#[derive(Default, Clone)]
pub struct ScoreDistribution {
    score_counts: Vec<usize>,
}

impl ScoreDistribution {
    pub fn add(&mut self, score: usize) {
        if self.score_counts.len() < score {
            self.score_counts.resize(score + 1, 0);
        }
        self.score_counts[score] += 1;
    }

    pub fn histogram(&self) {
        let max = self.score_counts.iter().max().unwrap();
        let maxw = max.to_string().len();
        let width = 40usize;
        let scale = max.div_ceil(width);
        let scale = scale.max(1);

        for (label, value) in self.score_counts.iter().enumerate() {
            let blocks = value / scale;
            println!(
                "{:>4}b | {value:>W$} | {}",
                label,
                "∎".repeat(blocks),
                W = maxw
            );
        }
    }

    pub fn min(&self) -> usize {
        self.score_counts.iter().position(|c| *c != 0).unwrap_or(0)
    }

    pub fn total(&self) -> usize {
        self.score_counts.iter().sum()
    }

    pub fn below(&self, score: usize) -> f32 {
        let cnt = self.score_counts.iter().take(score + 1).sum::<usize>() as f32;
        cnt / self.total() as f32
    }

    pub fn dump(&self, out: &mut impl Write) -> anyhow::Result<()> {
        for (sc, cnt) in self
            .score_counts
            .iter()
            .enumerate()
            .filter(|(_, c)| **c != 0)
        {
            writeln!(out, "{sc}: {cnt}")?;
        }

        Ok(())
    }
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let mut out: Box<dyn Write> = match args.out_path {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(stdout())),
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let prefilter_db = SplitDb::from_path(args.dir_path.join("prefilterDB"))
        .context("failed to build prefilter db")?;

    let query_map = Db::from_path(args.dir_path.join("queryDB_h"))
        .context("failed to build query db header")?
        .into_map();

    let mut queries = query_map.iter().collect::<Vec<_>>();
    queries.sort_by_key(|x| x.1);

    let queries = queries
        .into_iter()
        .map(|(q, _)| q.clone())
        .collect::<Vec<_>>();

    let target_map = Db::from_path(args.dir_path.join("targetDB_h"))
        .context("failed to build target db header")?
        .into_map();

    // ---

    let reader = BufReader::new(File::open(args.ga_hits)?);
    let mut sets = vec![HashSet::new(); query_map.len()];
    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            continue;
        }

        let tokens = line.split_whitespace().collect::<Vec<_>>();

        let query = tokens[0];
        let target = tokens[1];

        // let ga_sc = tokens[2].parse::<f32>()?;
        let nail_cloud_sc = tokens[4].parse::<f32>().ok();
        // let nail_sc = tokens[6].parse::<f32>().ok();
        // let hmmer_sc = tokens[10].parse::<f32>().ok();

        let q_idx = *query_map.get(query).unwrap();
        let t_idx = *target_map.get(target).unwrap();

        if nail_cloud_sc.is_some() {
            sets[q_idx].insert(t_idx);
        }
    }

    // ---

    let mut pf_distributions = vec![ScoreDistribution::default(); query_map.len()];
    let mut ga_distributions = vec![ScoreDistribution::default(); query_map.len()];

    #[allow(clippy::needless_range_loop)]
    for q_idx in 0..prefilter_db.len() {
        let pf_hits = prefilter_db.get(q_idx)?;

        let set = &sets[q_idx];
        let pf_dist = &mut pf_distributions[q_idx];
        let ga_dist = &mut ga_distributions[q_idx];
        for line in pf_hits.lines() {
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            let t_idx = tokens[0].parse::<usize>()?;
            let score = tokens[1].parse::<usize>()?;
            pf_dist.add(score);

            if set.contains(&t_idx) {
                ga_dist.add(score);
            }
        }

        writeln!(out, ">{}", queries[q_idx])?;
        pf_dist.dump(&mut out)?;
        writeln!(out, "//")?;
        ga_dist.dump(&mut out)?;
        writeln!(out, "//")?;
    }

    Ok(())
}
