use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::{env, path::Path};

use anyhow::{bail, Context};
use bioio::fasta::FastaRecord;
use bioio::{fasta::Fasta, stockholm::Stockholm};

use rand::seq::IndexedRandom;
use rand::{rng, Rng};
use rand_distr::{Distribution, Normal};

const N_BINS: usize = 41;
const BIN_START: usize = 10;
const BIN_END: usize = N_BINS + BIN_START - 1;

const MIN_LENGTH: usize = 200;

const DECOY_MEAN: f64 = 400.0;
const DECOY_STD_DEV: f64 = 50.0;
const DECOY_MIN_LEN: isize = 200;
const DECOY_MAX_LEN: isize = 600;
const DECOY_MIN_FRAG: usize = 2;
const DECOY_MAX_FRAG: usize = 4;

#[derive(Clone)]
struct Pair {
    /// The ID of a query/train MSA
    family: String,
    /// The ID of a sequence in the family MSA
    query: String,
    /// The name of a sequence in the target Fasta
    target: String,
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
        println!(
            "usage: pid <benchmark_dir> <query.sto> <target.sto> <source.sto> <source.fa> <n_sample> <decoy_ratio>"
        );
        return Ok(());
    };

    let bm_dir = Path::new(&args[1]);
    if !bm_dir.is_dir() {
        bail!("bm_dir: {bm_dir:?} is not a directory");
    }

    let mut query_sto = Stockholm::from_path(&args[2]).context("failed to parse query sto")?;
    let mut target_sto = Stockholm::from_path(&args[3]).context("failed to parse query sto")?;
    let src_sto = Stockholm::from_path(&args[4]).context("failed to parse source sto")?;
    let src_fa = Fasta::from_path(&args[5]).context("failed to parse source sto")?;
    let n_sample: usize = args[6].parse().context("failed to parse n_sample")?;
    let decoy_ratio: usize = args[7].parse().context("failed to parse decoy ratio")?;

    // what: filter all target/query seqs
    //       that are less than MIN_LENGTH
    //
    // why:  short sequences at very low %ID
    //       are basically impossible to find
    let len_filter = |sto: &mut Stockholm| {
        sto.records.values_mut().for_each(|rec| {
            let keep = rec
                .sequences
                .iter()
                .filter(|(_, seq)| seq.bytes().filter(|b| AMINO[*b as usize]).count() >= MIN_LENGTH)
                .map(|(name, _)| name.clone())
                .collect::<HashSet<_>>();

            rec.gs_meta.retain(|name, _| keep.contains(name));
            rec.sequences.retain(|name, _| keep.contains(name));
        });

        sto.records.retain(|_, rec| !rec.sequences.is_empty());
    };

    len_filter(&mut query_sto);
    len_filter(&mut target_sto);

    // what: remove families from the target set that don't
    //       appear in the query set (and vice versa)
    //
    // why:  if there's no available queries for a target,
    //       we'll never find it, and if there's no targets
    //       for a query, we're just wasting search time
    target_sto
        .records
        .retain(|fam, _| query_sto.records.contains_key(fam));

    query_sto
        .records
        .retain(|fam, _| target_sto.records.contains_key(fam));

    assert_eq!(
        target_sto.records.len(),
        query_sto.records.len(),
        "target/query family count mismatch"
    );

    // what: for each target seq, find the query seq
    //       with the highest pairwise %ID
    //
    // why:  we are building our query set such that
    //       each target gets its best possible query
    let mut pairs_by_bin: Vec<Vec<Pair>> = vec![vec![]; BIN_END + 1];
    for fam in query_sto.records.keys() {
        let src_seqs = &src_sto
            .get(fam)
            .map(|rec| &rec.sequences)
            .context("failed to retrieve seq family: \"{seq}\" from source stockholm")?;

        let query_seqs = src_seqs
            .iter()
            .filter(|(name, _)| {
                query_sto
                    .get(fam)
                    .expect("failed to retrieve seq family: \"{seq}\" from query stockholm")
                    .sequences
                    .contains_key(*name)
            })
            .collect::<Vec<_>>();

        let target_seqs = src_seqs
            .iter()
            .filter(|(name, _)| {
                target_sto
                    .get(fam)
                    .expect("failed to retrieve seq family: \"{seq}\" from target stockholm")
                    .sequences
                    .contains_key(*name)
            })
            .collect::<Vec<_>>();

        for (t_name, t_seq) in target_seqs.iter() {
            let mut best_pid = 0.0;
            let mut best_query = "";
            for (q_name, q_seq) in query_seqs.iter() {
                let pid = compute_pid(t_seq, q_seq);
                assert!(pid <= 0.5, "unexpected high identity: {}%", pid * 100.0);

                if pid > best_pid {
                    best_pid = pid;
                    best_query = q_name;
                }
            }
            let bin = (best_pid * 100.0).round() as usize;

            pairs_by_bin[bin].push(Pair {
                family: fam.clone(),
                query: best_query.to_string(),
                target: t_name.to_string(),
            })
        }
    }

    let mut target_writer = BufWriter::new(
        File::create(bm_dir.join("target.fa")).context("failed to open target.fa for output")?,
    );
    let mut tbl_writer = BufWriter::new(
        File::create(bm_dir.join("benchmark.tbl"))
            .context("failed to open benchmark.tbl for output")?,
    );
    writeln!(tbl_writer, "#identity family target query")?;

    // store queries in a hash to prevent duplicates, since it's
    // possible that two targets share a most similar query
    let mut queries: HashSet<FastaRecord> = HashSet::new();
    let mut rng = rng();

    println!("sampling {n_sample} from each bin...");
    #[allow(clippy::needless_range_loop)]
    for bin in BIN_START..=BIN_END {
        let bin_size = pairs_by_bin[bin].len();
        let pairs: Vec<&Pair> = if bin_size >= n_sample {
            pairs_by_bin[bin]
                .choose_multiple(&mut rng, n_sample)
                .collect()
        } else {
            println!("warning: bin {bin}% has {bin_size:>3}/{n_sample} samples",);
            pairs_by_bin[bin].iter().collect()
        };

        let extract_fn = |sto: &Stockholm, fam: &str, seq: &str| {
            sto.get(fam)
                .and_then(|r| r.get(seq))
                .map(|s| s.replace(['-', '.'], ""))
                .context("failed to extract seq from source stockholm")
        };

        for (pair_idx, pair) in pairs.iter().enumerate() {
            let query = extract_fn(&query_sto, &pair.family, &pair.query)
                .map(|seq| FastaRecord {
                    name: format!("{}|{}", pair.family, pair.query),
                    extra: "".to_string(),
                    seq,
                })
                .context("failed to build query fasta record")?;

            let target = extract_fn(&target_sto, &pair.family, &pair.target)
                .map(|seq| FastaRecord {
                    name: format!("{}|{}|{}%:{}", pair.family, pair.target, bin, pair_idx),
                    extra: "".to_string(),
                    seq,
                })
                .context("failed to build target fasta record")?;

            queries.insert(query);

            writeln!(target_writer, "{target}").context("failed to write to target.fa")?;

            writeln!(
                tbl_writer,
                "{}% {} {} {}",
                bin, pair.family, pair.target, pair.query
            )
            .context("failed to write to benchmark.tbl")?;
        }
    }

    let mut query_fa_writer = BufWriter::new(
        File::create(bm_dir.join("query.fa")).context("failed to open query.fa for output")?,
    );

    let mut queries = queries.into_iter().collect::<Vec<_>>();
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
        .collect::<HashSet<_>>();

    query_sto
        .records
        .iter()
        .filter(|(fam, _)| query_names.contains(&fam.as_str()))
        .try_for_each(|(_, r)| writeln!(query_sto_writer, "{r}"))?;

    let n_src = src_fa.records.len();
    let n_decoys = n_sample * N_BINS * decoy_ratio;

    println!("sampling {n_decoys} decoys from a source of {n_src} seqs...");

    let len_distribution = Normal::new(DECOY_MEAN, DECOY_STD_DEV).unwrap();
    for decoy_idx in 0..n_decoys {
        let decoy_len = (len_distribution.sample(&mut rng).round() as isize)
            .clamp(DECOY_MIN_LEN, DECOY_MAX_LEN) as usize;
        let frag_cnt: usize = rng.random_range(DECOY_MIN_FRAG..=DECOY_MAX_FRAG);
        let mut points: Vec<f64> = (0..frag_cnt).map(|_| rng.random_range(0.0..=1.0)).collect();

        points.append(&mut vec![0.0, 1.0]);
        points.sort_by(|a, b| a.partial_cmp(b).expect("NaN encountered"));

        let mut decoy_bytes = Vec::new();
        for frag_len in points
            .windows(2)
            .map(|p| ((p[1] - p[0]) * decoy_len as f64).round() as usize)
        {
            let mut n_to_add = frag_len;
            while n_to_add > 0 {
                let src_bytes = src_fa
                    .records
                    .get_index(rng.random_range(0..n_src))
                    .context("bad fasta index")?
                    .1
                    .seq
                    .as_bytes();

                let byte_sample = src_bytes.choose_multiple(&mut rng, n_to_add).copied();

                n_to_add -= byte_sample.len();
                decoy_bytes.extend(byte_sample);
            }
        }

        let decoy = FastaRecord {
            name: format!("decoy{decoy_idx}"),
            extra: "".to_string(),
            seq: std::str::from_utf8(&decoy_bytes)
                .context("failed to build decoy seq")?
                .to_string(),
        };

        writeln!(target_writer, "{decoy}").context("failed to write to target.fa")?;
    }

    Ok(())
}

//   -------------------------------------
//   distribution of lengths in Swissprot:
//   -------------------------------------
//   2 ..   22 [  3986 ]: ∎∎∎∎∎∎∎
//  22 ..   42 [  6032 ]: ∎∎∎∎∎∎∎∎∎∎
//  42 ..   62 [  9111 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
//  62 ..   82 [ 17284 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
//  82 ..  102 [ 23248 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 102 ..  122 [ 22480 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 122 ..  142 [ 24662 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 142 ..  162 [ 27948 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 162 ..  182 [ 21797 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 182 ..  202 [ 24136 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 202 ..  222 [ 25015 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 222 ..  242 [ 22346 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 242 ..  262 [ 22755 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 262 ..  282 [ 20434 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 282 ..  302 [ 20847 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 302 ..  322 [ 22050 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 322 ..  342 [ 20781 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 342 ..  362 [ 20953 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 362 ..  382 [ 18927 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 382 ..  402 [ 16276 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 402 ..  422 [ 14426 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 422 ..  442 [ 16088 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 442 ..  462 [ 14004 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 462 ..  482 [ 12909 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 482 ..  502 [ 11433 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 502 ..  522 [ 10788 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 522 ..  542 [  7193 ]: ∎∎∎∎∎∎∎∎∎∎∎∎
// 542 ..  562 [  7997 ]: ∎∎∎∎∎∎∎∎∎∎∎∎∎∎
// 562 ..  582 [  6409 ]: ∎∎∎∎∎∎∎∎∎∎∎
// 582 ..  602 [  5719 ]: ∎∎∎∎∎∎∎∎∎∎
// 602 ..  622 [  5505 ]: ∎∎∎∎∎∎∎∎∎
// 622 ..  642 [  5357 ]: ∎∎∎∎∎∎∎∎∎
// 642 ..  662 [  4196 ]: ∎∎∎∎∎∎∎
// 662 ..  682 [  3663 ]: ∎∎∎∎∎∎
// 682 ..  702 [  3836 ]: ∎∎∎∎∎∎
// 702 ..  722 [  3513 ]: ∎∎∎∎∎∎
// 722 ..  742 [  3030 ]: ∎∎∎∎∎
// 742 ..  762 [  2669 ]: ∎∎∎∎
// 762 ..  782 [  2214 ]: ∎∎∎
// 782 ..  802 [  2144 ]: ∎∎∎
// 802 ..  822 [  2081 ]: ∎∎∎
// 822 ..  842 [  1918 ]: ∎∎∎
// 842 ..  862 [  2069 ]: ∎∎∎
// 862 ..  882 [  2297 ]: ∎∎∎∎
// 882 ..  902 [  1930 ]: ∎∎∎
// 902 ..  922 [  1709 ]: ∎∎∎
// 922 ..  942 [  1494 ]: ∎∎
// 942 ..  962 [  1664 ]: ∎∎
// 962 ..  982 [  1263 ]: ∎∎
// 982 .. 1002 [   845 ]: ∎
