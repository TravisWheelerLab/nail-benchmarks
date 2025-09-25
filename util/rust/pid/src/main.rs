use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::{env, path::Path};

use anyhow::{bail, Context};
use bioio::fasta::FastaRecord;
use bioio::{fasta::Fasta, stockholm::Stockholm};

use rand::seq::{IndexedRandom, SliceRandom};
use rand::{rng, Rng};

const DECOY_RATIO: usize = 100;

const PID_MAX: usize = 25;

const SEQ_LEN_MIN: usize = 150;

const FAM_MIN: usize = 5;
const FAM_MAX: usize = 10;

#[derive(Clone)]
struct Pair {
    pid: usize,
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

    if args.len() < 6 {
        println!(
            "usage: pid <benchmark_dir> <query.sto> <target.sto> <source.sto> <source.fa> [max_pairs]"
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
    let afa_dir = bm_dir.join("afa/");

    let max_pairs: Option<usize> = if args.len() == 7 {
        Some(args[6].parse().context("failed to parse n_sample")?)
    } else {
        None
    };

    assert_eq!(
        target_sto.records.len(),
        query_sto.records.len(),
        "target/query family count mismatch"
    );

    println!("{} sequence families found", target_sto.records.len());

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
                .filter(|(_, seq)| {
                    seq.bytes().filter(|b| AMINO[*b as usize]).count() >= SEQ_LEN_MIN
                })
                .map(|(name, _)| name.clone())
                .collect::<HashSet<_>>();

            rec.gs_meta.retain(|name, _| keep.contains(name));
            rec.sequences.retain(|name, _| keep.contains(name));
        });

        sto.records.retain(|_, rec| !rec.sequences.is_empty());
    };

    len_filter(&mut query_sto);
    len_filter(&mut target_sto);

    // only keep query fams that still have at least 10 sequences
    query_sto.records.retain(|_, rec| rec.sequences.len() >= 10);

    // what: remove families from the target set that don't
    //       appear in the query set (and vice versa)
    //
    // why:  if there's no available queries for a target,
    //       we'll never find it, and if there's no targets
    //       for a query, we're not looking for anything
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

    println!(
        "{} sequence families remain after length filter (>={SEQ_LEN_MIN})",
        target_sto.records.len()
    );

    // what: for each target seq, find the query seq
    //       with the highest pairwise %ID
    //
    // why:  we are building our query set such that
    //       each target gets its best possible query
    let mut pairs_by_fam: HashMap<String, Vec<Pair>> = query_sto
        .records
        .keys()
        .map(|fam| (fam.clone(), vec![]))
        .collect();

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

        let mut keep = HashSet::new();
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

            if bin <= PID_MAX {
                keep.insert(t_name);
                pairs_by_fam
                    .get_mut(fam)
                    .context("no pair vec for fam")?
                    .push(Pair {
                        pid: bin,
                        family: fam.clone(),
                        query: best_query.to_string(),
                        target: t_name.to_string(),
                    })
            }
        }

        let sto = target_sto.get_mut(fam).context("")?;

        sto.sequences.retain(|n, _| keep.contains(&n));
        sto.gs_meta.retain(|n, _| keep.contains(&n));
    }

    pairs_by_fam.retain(|_, pairs| pairs.len() > FAM_MIN);

    pairs_by_fam.iter_mut().for_each(|(_, pairs)| {
        pairs.sort_by(|a, b| a.pid.cmp(&b.pid));
        pairs.truncate(FAM_MAX);
    });

    let mut rng = rng();

    let mut pairs: Vec<Pair> = pairs_by_fam.into_values().flatten().collect();
    pairs = match max_pairs {
        Some(max) => pairs.choose_multiple(&mut rng, max).cloned().collect(),
        None => pairs,
    };
    pairs.sort_by(|a, b| a.pid.cmp(&b.pid));

    let mut tbl_writer = BufWriter::new(
        File::create(bm_dir.join("benchmark.tbl"))
            .context("failed to open benchmark.tbl for output")?,
    );
    writeln!(tbl_writer, "#identity family target query")?;

    let mut targets: Vec<FastaRecord> = Vec::new();
    // store queries in a hash to prevent duplicates, since it's
    // possible that two targets share a most similar query
    let mut queries: HashSet<FastaRecord> = HashSet::new();

    let extract_fn = |sto: &Stockholm, fam: &str, seq: &str| {
        sto.get(fam)
            .and_then(|r| r.get(seq))
            .map(|s| s.replace(['-', '.'], ""))
            .context("failed to extract seq from stockholm")
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
                name: format!("{}|{}|{}%:{}", pair.family, pair.target, pair.pid, pair_idx),
                extra: "".to_string(),
                seq,
            })
            .context("failed to build target fasta record")?;

        targets.push(target);
        queries.insert(query);

        writeln!(
            tbl_writer,
            "{}% {} {} {}",
            pair.pid, pair.family, pair.target, pair.query
        )
        .context("failed to write to benchmark.tbl")?;
    }

    println!("done");

    let mut target_writer = BufWriter::new(
        File::create(bm_dir.join("target.fa")).context("failed to open target.fa for output")?,
    );

    targets
        .iter()
        .try_for_each(|t| writeln!(target_writer, "{t}").context("failed to write to target.fa"))?;

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

    fs::create_dir(&afa_dir)?;

    let query_names = queries
        .iter()
        .map(|q| q.name.split('|').next().unwrap())
        .collect::<HashSet<_>>();

    query_sto
        .records
        .retain(|fam, _| query_names.contains(&fam.as_str()));

    query_sto.records.iter().try_for_each(|(fam, rec)| {
        writeln!(query_sto_writer, "{rec}")?;
        let mut afa_writer = BufWriter::new(File::create(afa_dir.join(format!("{fam}.afa")))?);
        rec.sequences.iter().try_for_each(|(name, seq)| {
            writeln!(afa_writer, ">{fam}:{name}")?;
            writeln!(afa_writer, "{seq}")
        })
    })?;

    let n_src = src_fa.records.len();
    let n_decoys = pairs.len() * DECOY_RATIO;

    println!("sampling {n_decoys} decoys from a source of {n_src} seqs...");

    let mut src_bytes: &[u8] = &[];
    for decoy_idx in 0..n_decoys {
        // pick a length by taking a random target
        let decoy_len = targets
            .get(rng.random_range(0..targets.len()))
            .context("bad fasta index")?
            .seq
            .len();

        // grab from the source seqs until we
        // find one of at least that length
        while src_bytes.len() < decoy_len {
            src_bytes = src_fa
                .records
                .get_index(rng.random_range(0..n_src))
                .context("bad fasta index")?
                .1
                .seq
                .as_bytes();
        }

        let max_start = src_bytes.len() - decoy_len;
        let start = rng.random_range(0..=max_start);
        let end = start + decoy_len;

        let mut sample: Vec<u8> = src_bytes[start..end].to_vec();
        sample.shuffle(&mut rng);

        let decoy = FastaRecord {
            name: format!("decoy{decoy_idx}"),
            extra: "".to_string(),
            seq: std::str::from_utf8(&sample)
                .context("failed to build decoy seq")?
                .to_string(),
        };

        writeln!(target_writer, "{decoy}").context("failed to write to target.fa")?;

        src_bytes = &[]
    }
    println!("done");

    Ok(())
}
