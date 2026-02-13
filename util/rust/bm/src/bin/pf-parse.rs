use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::{Arc, Mutex},
};

use anyhow::Context;
use bioio::{
    mmseqs::Db,
    tbl::{HitTable, NailTable},
};
use clap::{Args, Parser, Subcommand};
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    ThreadPoolBuilder,
};

#[derive(Subcommand)]
enum SubCommands {
    Prog(ProgArgs),
    Score(ScoreArgs),
}

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    pub command: SubCommands,
}

#[derive(Args)]
struct CommonArgs {
    #[arg(value_name = "prefilterDB")]
    pdb_path: PathBuf,

    #[arg(value_name = "alignDB")]
    adb_path: PathBuf,

    #[arg(value_name = "queryDB_h")]
    qdb_h_path: PathBuf,

    #[arg(value_name = "targetDB_h")]
    tdb_h_path: PathBuf,

    // #[arg(value_name = "query.hmm")]
    // query_path: PathBuf,
    #[arg(value_name = "results.tbl")]
    nail_tbl_path: PathBuf,

    #[arg(value_name = "path")]
    out_path: PathBuf,

    #[arg(short = 't', default_value_t = 4usize, value_name = "N")]
    num_threads: usize,
}

#[derive(Args)]
struct ProgArgs {
    #[command(flatten)]
    common: CommonArgs,

    #[arg(long, value_name = "n", default_value_t = 200usize)]
    prog_n: usize,
}

#[derive(Args)]
struct ScoreArgs {
    #[command(flatten)]
    common: CommonArgs,
}

fn set_threads(num_threads: usize) -> anyhow::Result<()> {
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("failed to build rayon global threadpool")
}

fn main() -> anyhow::Result<()> {
    let now = std::time::Instant::now();
    match Cli::parse().command {
        SubCommands::Prog(args) => prog(args)?,
        SubCommands::Score(args) => score(args)?,
    };
    println!("took: {:?}", now.elapsed());
    Ok(())
}

struct Hit {
    psc: usize,
    asc: Option<usize>,
    nsc: Option<f64>,
}

fn hits(
    pdb_records: String,
    adb_records: String,
    nail_records: Option<&HashMap<String, f64>>,
    target_names: &[String],
) -> anyhow::Result<Vec<Hit>> {
    let mut map = HashMap::with_capacity(pdb_records.len() / 10);

    for line in pdb_records.lines() {
        let mut tokens = line.split_whitespace();
        let tid: usize = tokens
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing tid"))?
            .parse()?;
        let psc: usize = tokens
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing psc"))?
            .parse()?;

        map.insert(
            tid,
            Hit {
                psc,
                asc: None,
                nsc: None,
            },
        );
    }

    for line in adb_records.lines() {
        let mut tokens = line.split_whitespace();
        let tid: usize = tokens
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing tid"))?
            .parse()?;
        let asc: usize = tokens
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing asc"))?
            .parse()?;

        map.get_mut(&tid).expect("no tid for align record").asc = Some(asc);
    }

    if let Some(recs) = nail_records {
        map.iter_mut().for_each(|(&tid, hit)| {
            let tname = &target_names[tid];
            hit.nsc = recs.get(tname).cloned();
        });
    }

    Ok(map.into_values().collect())
}

fn score(args: ScoreArgs) -> anyhow::Result<()> {
    set_threads(args.common.num_threads)?;

    if let Some(parent) = args.common.out_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let out = Arc::new(Mutex::new(BufWriter::new(File::create(
        args.common.out_path,
    )?)));

    let mut pdb = Db::from_path(args.common.pdb_path).context("failed to build prefilter db")?;
    let mut adb = Db::from_path(args.common.adb_path).context("failed to build adb")?;
    let mut qdb_h = Db::from_path(args.common.qdb_h_path).context("failed to build qdb_h")?;
    let mut tdb_h = Db::from_path(args.common.tdb_h_path).context("failed to build qdb_h")?;

    let mut target_names: Vec<String> = vec![];

    for idx in 0..tdb_h.len() {
        target_names.push(tdb_h.get(idx)?.trim_end().to_string());
    }

    let tbl = HitTable::from_path::<_, NailTable>(args.common.nail_tbl_path)?;
    let mut nail_map: HashMap<String, HashMap<String, f64>> = HashMap::new();

    tbl.hits.into_iter().for_each(|h| {
        nail_map
            .entry(h.query)
            .or_default()
            .insert(h.target, h.score);
    });

    (0..qdb_h.len())
        .map(|idx| {
            let mut q = qdb_h.get(idx).unwrap();
            q = q.trim().to_string();
            let p = pdb.get(idx).unwrap();
            let a = adb.get(idx).unwrap();
            let n = nail_map.get(&q);
            (q, p, a, n, out.clone())
        })
        .par_bridge()
        .for_each(|(query, pdb_records, adb_records, nail_records, out)| {
            let mut hits = hits(pdb_records, adb_records, nail_records, &target_names).unwrap();

            hits.sort_by(|a, b| b.psc.cmp(&a.psc));

            let mut buf = vec![];

            if !hits.is_empty() {
                write!(buf, "{query}").unwrap();

                for chunk in hits.chunk_by(|a, b| a.psc == b.psc) {
                    let psc = chunk[0].psc;

                    let hit_cnt = chunk.len();
                    let seed_cnt = chunk.iter().filter(|h| h.asc.is_some()).count();
                    let seed_frac = seed_cnt as f32 / chunk.len() as f32;

                    let nail_cnt = chunk.iter().filter(|h| h.nsc.is_some()).count();
                    let nail_frac = nail_cnt as f32 / chunk.len() as f32;

                    write!(
                        buf,
                        ",({},{:.4},{},{:.4},{},{})",
                        psc, seed_frac, seed_cnt, nail_frac, nail_cnt, hit_cnt
                    )
                    .unwrap();
                }

                writeln!(buf).unwrap();
            }

            match out.lock() {
                Ok(mut guard) => {
                    guard.write_all(&buf).unwrap();
                    guard.flush().unwrap();
                }
                Err(_) => panic!(),
            }
        });

    Ok(())
}

fn prog(args: ProgArgs) -> anyhow::Result<()> {
    // if let Some(parent) = args.common.out_path.parent() {
    //     std::fs::create_dir_all(parent)?;
    // }
    // let mut out = BufWriter::new(File::create(args.common.out_path)?);
    // let mut pdb = Db::from_path(args.common.pdb_path).context("failed to build prefilter db")?;
    // let mut qdb_h = Db::from_path(args.common.qdb_h_path).context("failed to build qdb_h")?;
    // let mut adb = Db::from_path(args.common.adb_path).context("failed to build adb")?;

    // for idx in 0..qdb_h.len() {
    //     let mut query = qdb_h.get(idx)?;
    //     query = query.trim().to_string();

    //     let mut hits = hits(idx, &mut pdb, &mut adb)?;

    //     hits.sort_by(|a, b| b.psc.cmp(&a.psc));

    //     if !hits.is_empty() {
    //         write!(out, "{query}")?;

    //         let mut n = args.prog_n;
    //         let g = 1usize;
    //         let mut tot_hit_cnt = 0;
    //         let mut tot_seed_cnt = 0;
    //         let mut bin_start = 0;
    //         let mut iter = hits.iter().peekable();
    //         while iter.peek().is_some() {
    //             let chunk = iter.by_ref().take(n).collect::<Vec<_>>();
    //             let bin_seed_cnt = chunk.iter().filter(|h| h.asc.is_some()).count();

    //             tot_hit_cnt += chunk.len();
    //             tot_seed_cnt += bin_seed_cnt;

    //             let asc_min = chunk.iter().map(|h| h.psc).min().unwrap();
    //             let asc_max = chunk.iter().map(|h| h.psc).max().unwrap();
    //             let asc_avg =
    //                 chunk.iter().map(|h| h.psc).sum::<usize>() as f32 / chunk.len() as f32;

    //             let bin_end = bin_start + chunk.len();

    //             let bin_seed_frac = bin_seed_cnt as f32 / chunk.len() as f32;
    //             let tot_seed_frac = tot_seed_cnt as f32 / tot_hit_cnt as f32;

    //             write!(
    //                 out,
    //                 ",({},{},{},{},{:.1},{:.4},{:.4})",
    //                 bin_start, bin_end, asc_min, asc_max, asc_avg, tot_seed_frac, bin_seed_frac
    //             )?;

    //             bin_start = bin_end;
    //             n *= g;
    //         }

    //         writeln!(out)?;
    //     }
    // }

    Ok(())
}
