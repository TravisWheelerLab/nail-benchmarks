use std::{
    collections::HashMap,
    fs::File,
    io::{stdout, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use clap::{Parser, Subcommand};
use glob::glob;
use rayon::prelude::*;

// | GA |    GA     |  nail  |    nail   |  nail  |   nail    | mmseqs |  mmseqs   | hmmer |  hmmer    | dom | dom | sig |
// | sc |  P-value  | cld sc |cld P-value|   sc   |  P-value  |   sc   |  P-value  |  sc   |  P-value  | sum | max | dom | dom scores
//  ---- ----------- -------- ----------- -------- ----------- -------- ----------- ------- ----------- ----- ----- ----- ------------

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    Dist {
        #[arg(value_name = "hits/")]
        dir_path: PathBuf,

        #[arg(short, long, value_name = "path")]
        out_path: Option<PathBuf>,

        #[arg(short, value_name = "n", default_value_t = 8)]
        threads: usize,
    },
    Ga {
        #[arg(value_name = "hits.tbl")]
        tbl_path: PathBuf,
    },
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().cmd {
        Cmd::Dist {
            dir_path,
            out_path,
            threads,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            let paths = glob(dir_path.join("*.tbl").to_str().context("invalid glob")?)?
                .filter_map(Result::ok)
                .collect::<Vec<_>>();

            let mut out: Box<dyn Write> = match out_path {
                Some(p) => Box::new(BufWriter::new(File::create(p)?)),
                None => Box::new(BufWriter::new(stdout())),
            };

            query_distribution(paths, &mut out)?;
        }
        Cmd::Ga { tbl_path } => {
            max_seqs_info(tbl_path)?;
        }
    }

    Ok(())
}

fn query_distribution(paths: Vec<PathBuf>, out: &mut dyn Write) -> anyhow::Result<()> {
    let maps = paths
        .par_iter()
        .map(|path| {
            let mut map: HashMap<String, usize> = HashMap::new();
            let reader = BufReader::new(File::open(path).unwrap());

            for line in reader.lines() {
                let line = line.unwrap();

                if line.starts_with('#') {
                    continue;
                }

                let tokens = line.split_whitespace().collect::<Vec<_>>();

                let query = tokens[0];
                let ga_sc = tokens[2].parse::<f32>().unwrap();
                let hmmer_sc = tokens[10].parse::<f32>().ok();

                if let Some(sc) = hmmer_sc {
                    if sc >= ga_sc {
                        match map.get_mut(query) {
                            Some(cnt) => *cnt += 1,
                            None => {
                                map.insert(query.to_string(), 1);
                            }
                        }
                    }
                }
            }

            map
        })
        .collect::<Vec<_>>();

    let mut dist: Vec<(String, Vec<usize>)> = vec![];

    fn cv(counts: &[usize]) -> f32 {
        let n = counts.len() as f32;
        let mean = counts.iter().map(|&x| x as f32).sum::<f32>() / n;
        let var = counts
            .iter()
            .map(|&x| {
                let d = x as f32 - mean;
                d * d
            })
            .sum::<f32>()
            / n;
        var.sqrt() / (mean)
    }

    for k in maps[0].keys().cloned() {
        let mut v: Vec<usize> = vec![];
        for m in &maps {
            match m.get(&k) {
                Some(cnt) => v.push(*cnt),
                None => v.push(0),
            }
        }
        dist.push((k, v));
    }

    dist.sort_by_key(|k| *k.1.iter().max().unwrap());

    for (q, v) in dist.iter() {
        writeln!(out, "{q:20} {:.3} {v:?}", cv(v))?
    }

    Ok(())
}

#[derive(Default)]
struct MaxSeqsRecord {
    hmmer_hits: usize,
    hmmer_ga_hits: usize,
    nail_seeds: usize,
    nail_hits: usize,
    nail_ga_hits: usize,
}

fn max_seqs_info(path: impl AsRef<Path>) -> anyhow::Result<()> {
    let reader = BufReader::new(File::open(path)?);

    let mut map: HashMap<String, MaxSeqsRecord> = HashMap::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            continue;
        }

        let tokens = line.split_whitespace().collect::<Vec<_>>();

        let query = tokens[0];
        let ga_sc = tokens[2].parse::<f32>()?;
        let nail_cloud_sc = tokens[4].parse::<f32>().ok();
        let nail_sc = tokens[6].parse::<f32>().ok();
        let hmmer_sc = tokens[10].parse::<f32>().ok();

        let record = match map.get_mut(query) {
            Some(r) => r,
            None => {
                map.insert(query.to_string(), MaxSeqsRecord::default());
                map.get_mut(query).unwrap()
            }
        };

        if let Some(sc) = hmmer_sc {
            record.hmmer_hits += 1;
            record.hmmer_ga_hits += (sc >= ga_sc) as usize
        }

        if nail_cloud_sc.is_some() {
            record.nail_seeds += 1;
        }

        if let Some(sc) = nail_sc {
            record.nail_hits += 1;
            record.nail_ga_hits += (sc >= ga_sc) as usize
        }
    }

    let mut records = map.into_iter().collect::<Vec<_>>();
    records.sort_by(|a, b| b.1.hmmer_ga_hits.cmp(&a.1.hmmer_ga_hits));

    let mut w1 = 0;
    let mut w2 = 0;
    let mut w3 = 0;
    let mut w4 = 0;
    let mut w5 = 0;
    let mut w6 = 0;

    records.iter().for_each(|(q, r)| {
        w1 = w1.max(q.len());
        w2 = w2.max(r.hmmer_hits.to_string().len());
        w3 = w3.max(r.hmmer_ga_hits.to_string().len());
        w4 = w4.max(r.nail_seeds.to_string().len());
        w5 = w5.max(r.nail_hits.to_string().len());
        w6 = w5.max(r.nail_ga_hits.to_string().len());
    });

    for (query, record) in records {
        println!(
            "{query:W1$} {:W2$} {:W3$} {:W4$} {:W5$} {:W6$} ",
            record.hmmer_hits,
            record.hmmer_ga_hits,
            record.nail_seeds,
            record.nail_hits,
            record.nail_ga_hits,
            W1 = w1,
            W2 = w2,
            W3 = w3,
            W4 = w4,
            W5 = w5,
            W6 = w6,
        );
    }

    Ok(())
}
