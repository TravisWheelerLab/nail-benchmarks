use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::{collections::HashMap, env, path::Path};

use anyhow::{bail, Context};
use bioio::fasta::FastaRecord;
use bioio::{fasta::Fasta, stockholm::Stockholm};

use rand::rng;
use rand::seq::index::sample;

const N_BINS: usize = 41;
const BIN_START: usize = 10;
const BIN_MAX: usize = N_BINS + BIN_START - 1;

const MIN_LENGTH: usize = 200;

#[derive(Clone)]
struct Pair {
    /// The ID of a query/train MSA
    family: String,
    /// The ID of a sequence in the family MSA
    query: String,
    /// The name of a sequence in the target Fasta
    target: String,
}

struct Target {
    /// The name of the Fasta record
    name: String,
    /// The name of the domain planted in the sequence
    domain: String,
}

impl Display for Pair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} - {}", self.target, self.query)
    }
}

static AMINO: [bool; 256] = {
    let mut t = [false; 256];
    t[b'A' as usize] = true;
    t[b'C' as usize] = true;
    t[b'D' as usize] = true;
    t[b'E' as usize] = true;
    t[b'F' as usize] = true;
    t[b'G' as usize] = true;
    t[b'H' as usize] = true;
    t[b'I' as usize] = true;
    t[b'K' as usize] = true;
    t[b'L' as usize] = true;
    t[b'M' as usize] = true;
    t[b'N' as usize] = true;
    t[b'P' as usize] = true;
    t[b'Q' as usize] = true;
    t[b'R' as usize] = true;
    t[b'S' as usize] = true;
    t[b'T' as usize] = true;
    t[b'V' as usize] = true;
    t[b'W' as usize] = true;
    t[b'Y' as usize] = true;
    t[b'a' as usize] = true;
    t[b'c' as usize] = true;
    t[b'd' as usize] = true;
    t[b'e' as usize] = true;
    t[b'f' as usize] = true;
    t[b'g' as usize] = true;
    t[b'h' as usize] = true;
    t[b'i' as usize] = true;
    t[b'k' as usize] = true;
    t[b'l' as usize] = true;
    t[b'm' as usize] = true;
    t[b'n' as usize] = true;
    t[b'p' as usize] = true;
    t[b'q' as usize] = true;
    t[b'r' as usize] = true;
    t[b's' as usize] = true;
    t[b't' as usize] = true;
    t[b'v' as usize] = true;
    t[b'w' as usize] = true;
    t[b'y' as usize] = true;
    t
};

fn compute_pid(s1: &str, s2: &str) -> f32 {
    assert_eq!(s1.len(), s2.len());

    let mut match_cnt = 0;
    let mut pos_cnt = 0;

    s1.as_bytes()
        .iter()
        .zip(s2.as_bytes().iter())
        .for_each(|(&a, &b)| {
            if AMINO[a as usize] || AMINO[b as usize] {
                pos_cnt += 1;
                match_cnt += (a == b) as usize
            }
        });

    match_cnt as f32 / pos_cnt as f32
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 8 {
        println!("usage: pid <benchmark_dir> <positives.fa> <decoys.fa> <query.sto> <source.sto> <n_sample> <decoy_ratio>");
        return Ok(());
    };

    let bm_dir = Path::new(&args[1]);
    if !bm_dir.is_dir() {
        bail!("bm_dir: {bm_dir:?} is not a directory");
    }
    let positive_fa = Fasta::from_path(&args[2]).context("failed to parse positive fa")?;
    let decoy_fa = Fasta::from_path(&args[3]).context("failed to parse deocy fa")?;
    let query_sto = Stockholm::from_path(&args[4]).context("failed to parse query sto")?;
    let src_sto = Stockholm::from_path(&args[5]).context("failed to parse source sto")?;
    let n_sample: usize = args[6].parse().context("failed to parse n_sample")?;
    let decoy_ratio: usize = args[7].parse().context("failed to parse decoy ratio")?;

    let positive_fa: HashMap<String, FastaRecord> = positive_fa
        .records
        .into_iter()
        .filter(|(name, _)| {
            let mut range = name.split('/').next_back().expect("no range").split('-');
            let s: usize = range
                .next()
                .expect("no start")
                .parse()
                .expect("failed to parse start");
            let e: usize = range
                .next()
                .expect("no end")
                .parse()
                .expect("failed to parse end");

            e - s >= MIN_LENGTH
        })
        .collect();

    let mut targets_by_fam: HashMap<String, Vec<Target>> = HashMap::new();

    positive_fa.iter().for_each(|(name, rec)| {
        let name = name.clone();
        let fam = name
            .split('/')
            .next()
            .expect("failed to produce query family name")
            .to_string();

        let domain = rec
            .extra
            .split_whitespace()
            .last()
            .expect("failed to produce domain name")
            .to_string();

        targets_by_fam
            .entry(fam)
            .or_default()
            .push(Target { name, domain });
    });

    let mut pairs_by_bin: Vec<Vec<Pair>> = vec![vec![]; BIN_MAX + 1];

    targets_by_fam.iter().for_each(|(fam, targets)| {
        let query_msa = query_sto.records.get(fam).unwrap();
        let src_msa = src_sto.records.get(fam).unwrap();

        let src_names_in_query = query_msa
            .sequences
            .keys()
            .map(|s| s.to_owned())
            .collect::<Vec<_>>();

        let src_target_seqs: Vec<&String> = targets
            .iter()
            .map(|t| src_msa.sequences.get(&t.domain).unwrap())
            .collect();

        let src_query_seqs: Vec<&String> = src_names_in_query
            .iter()
            .map(|n| src_msa.sequences.get(n).unwrap())
            .collect();

        for (t_seq, target) in src_target_seqs
            .iter()
            .filter(|s| {
                s.as_bytes()
                    .iter()
                    .map(|b| *b as usize)
                    .filter(|b| AMINO[*b])
                    .sum::<usize>()
                    >= MIN_LENGTH
            })
            .zip(targets.iter())
        {
            let mut best_pid = 0.0;
            let mut best_q_name = "".to_string();

            for (q_seq, src_q_name) in src_query_seqs.iter().zip(src_names_in_query.iter()) {
                let pid = compute_pid(t_seq, q_seq);
                if pid > best_pid {
                    best_pid = pid;
                    best_q_name = src_q_name.to_string();
                }
            }

            let bin = (best_pid * 100.0).round() as usize;

            if bin > 50 {
                panic!("unexpected high %ID: {best_pid}");
            }

            pairs_by_bin[bin].push(Pair {
                family: fam.clone(),
                query: best_q_name,
                target: target.name.to_string(),
            })
        }
    });

    let mut target_writer = BufWriter::new(
        File::create(bm_dir.join("target.fa")).context("failed to open target.fa for output")?,
    );
    let mut tbl_writer = BufWriter::new(
        File::create(bm_dir.join("benchmark.tbl"))
            .context("failed to open benchmark.tbl for output")?,
    );
    writeln!(tbl_writer, "#target domain query")?;

    let mut queries: HashMap<String, FastaRecord> = HashMap::new();
    let mut rng = rng();

    // let min_bin_size = pairs_by_bin
    //     .iter()
    //     .skip(BIN_START)
    //     .map(|b| b.len())
    //     .min()
    //     .unwrap();

    // let n_sample = n_sample.min(min_bin_size).max(MIN_SAMPLE);
    println!("sampling {n_sample} from each bin...");

    (BIN_START..BIN_START + N_BINS)
        .try_for_each(|bin| {
            let length = pairs_by_bin[bin].len();

            let sample_indices = if length >= n_sample {
                let mut sample = sample(&mut rng, pairs_by_bin[bin].len(), n_sample).into_vec();
                sample.sort();
                sample
            } else {
                let l = pairs_by_bin[bin].len();
                println!("warning: bin {bin}% has {l:>3}/{n_sample} samples");
                (0..l).collect()
            };

            let pair_samples: Vec<&Pair> = sample_indices
                .into_iter()
                .map(|i| &pairs_by_bin[bin][i])
                .collect();

            pair_samples
                .iter()
                .enumerate()
                .try_for_each(|(pair_idx, p)| {
                    let query = FastaRecord {
                        name: format!("{}|{}", p.family, p.query),
                        extra: "".to_string(),
                        sequence: query_sto
                            .records
                            .get(&p.family)
                            .expect("")
                            .sequences
                            .get(&p.query)
                            .unwrap()
                            .replace(['-', '.'], ""),
                    };
                    match queries.get(&query.name) {
                        Some(r) => {
                            assert!(*r == query);
                        }
                        None => {
                            queries.insert(query.name.clone(), query);
                        }
                    };

                    let mut target = positive_fa.get(&p.target).unwrap().clone();
                    let domain_coords = target.name.split('/').nth(2).unwrap();
                    let domain = target.extra.split_whitespace().nth(1).unwrap();
                    target.name = format!("{}|{}|{}%:{}", p.family, domain_coords, bin, pair_idx);

                    writeln!(target_writer, "{target}").context("failed to write to target.fa")?;
                    writeln!(tbl_writer, "{} {} {}", target.name, domain, p.query)
                        .context("failed to write to benchmark.tbl")
                })
        })
        .context("bin sampling failed")?;

    if decoy_fa.records.is_empty() {
        println!("no decoys available");
    } else {
        let n_decoys = n_sample * N_BINS * decoy_ratio;
        println!("taking {n_decoys} decoys...");
        decoy_fa
            .records
            .values()
            .take(n_decoys)
            .try_for_each(|r| writeln!(target_writer, "{r}"))?;
    }

    let mut query_fa_writer = BufWriter::new(
        File::create(bm_dir.join("query.fa")).context("failed to open query.fa for output")?,
    );

    let mut queries: Vec<FastaRecord> = queries.into_values().collect();
    queries.sort_by(|a, b| a.name.cmp(&b.name));

    queries.iter().try_for_each(|q| {
        writeln!(query_fa_writer, "{q}").context("failed to write to query.fa")
    })?;

    let mut query_sto_writer = BufWriter::new(
        File::create(bm_dir.join("query.sto")).context("failed to open query.sto for output")?,
    );

    let query_names = queries
        .iter()
        .map(|q| q.name.split('|').next().unwrap())
        .collect::<Vec<_>>();

    query_sto
        .records
        .iter()
        .filter(|(fam, _)| query_names.contains(&fam.as_str()))
        .try_for_each(|(_, r)| writeln!(query_sto_writer, "{r}"))?;

    Ok(())
}
