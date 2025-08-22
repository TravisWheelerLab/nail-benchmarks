use std::{
    collections::HashMap,
    env,
    fs::File,
    path::{Path, PathBuf},
};

use anyhow::{anyhow, Context};
use bioio::{
    fasta::Fasta,
    tools::{hmmer, mmseqs, nail, Hit},
};

const N_BINS: usize = 41;
const BIN_START: usize = 10;
const FPR: f32 = 0.01;

fn extract_target_bin(target_name: &str) -> anyhow::Result<usize> {
    target_name
        .split("|")
        .nth(1)
        .ok_or(anyhow!("name doesn't have bin"))?
        .strip_suffix('%')
        .ok_or(anyhow!("no % suffix"))?
        .parse()
        .context("failed to parse bin")
}

struct TargetInfo<'a> {
    query: &'a str,
    start: usize,
    end: usize,
}

fn extract_target_info(hit: &Hit) -> anyhow::Result<TargetInfo> {
    let tokens = hit
        .target
        .split('|')
        .next()
        .ok_or(anyhow!("target name doesn't have a query"))?
        .split(':')
        .collect::<Vec<_>>();

    let range = tokens[1].split('-').collect::<Vec<_>>();
    Ok(TargetInfo {
        query: tokens[1],
        start: range[0].parse()?,
        end: range[1].parse()?,
    })
}

fn extract_query_name(hit: &Hit) -> anyhow::Result<&str> {
    let name = hit
        .query
        .split('|')
        .next()
        .ok_or(anyhow!("query doesn't have a name"))?;

    let name = match name.strip_suffix("-consensus") {
        Some(n) => n,
        None => name,
    };

    Ok(name)
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("usage: pid <benchmark_dir> [figures]");
        return Ok(());
    };

    let bm_dir = Path::new(&args[1]);
    let results_dir = bm_dir.join("results");

    let fig_dir = if args.len() > 2 {
        Path::new(&args[2])
    } else {
        Path::new("./figures/")
    };

    let target_fa = Fasta::from_path(bm_dir.join("target.fa"))?;
    let query_fa = Fasta::from_path(bm_dir.join("query.fa"))?;
    let cons_fa = Fasta::from_path(bm_dir.join("query.cons.fa"))?;

    let mut num_targets_by_bin = vec![0; BIN_START + N_BINS];
    query_fa
        .records
        .keys()
        .filter(|n| n.contains("decoy"))
        .try_for_each(|n| -> anyhow::Result<()> {
            let bin: usize = extract_target_bin(n)?;
            num_targets_by_bin[bin] += 1;
            Ok(())
        })?;

    let num_seq_queries = query_fa.records.len();
    let num_prf_queries = cons_fa.records.len();

    let seq_decoy_cnt = (FPR * num_seq_queries as f32) as usize;
    let prf_decoy_cnt = (FPR * num_prf_queries as f32) as usize;

    let mut hits = HashMap::new();
    let name_fn = |p: PathBuf| Some((p.clone(), p.file_stem()?.to_str()?.to_string()));

    for (path, name) in glob::glob(results_dir.join("hmmer*.domtbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.insert(name, hmmer::parse_domtbl(File::open(path)?)?);
    }

    for (path, name) in glob::glob(results_dir.join("mmseqs*.tbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.insert(name, mmseqs::parse_tbl(File::open(path)?)?);
    }

    for (path, name) in glob::glob(results_dir.join("nail*.tbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.insert(name, nail::parse_tbl(File::open(path)?)?);
    }

    for (name, list) in hits.iter_mut() {
        let n_decoys = match name.contains("seq") {
            true => seq_decoy_cnt,
            false => prf_decoy_cnt,
        };

        list.sort_by(|a, b| a.e_value.partial_cmp(&b.e_value).unwrap());

        let mut num_hits_by_bin = vec![0; BIN_START + N_BINS];

        let mut decoys_found = 0;
        for hit in list.iter() {
            match hit.target.contains("decoy") {
                true => {
                    decoys_found += 1;
                    if decoys_found >= n_decoys {
                        break;
                    }
                }
                false => {
                    let info = extract_target_info(hit)?;
                    let bin = extract_target_bin(&hit.target)?;
                    num_hits_by_bin[bin] += 1;
                }
            }
        }

        println!("{name} | {n_decoys} decoys");
        num_hits_by_bin
            .iter()
            .enumerate()
            .skip(BIN_START)
            .for_each(|(b, c)| println!("{b} | {c}"));
        println!();
    }

    Ok(())
}
