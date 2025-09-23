use std::{
    collections::{HashMap, HashSet},
    env,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use bioio::tbl::{BlastTable, Hit, HitTable, HmmerTable, NailTable};

use anyhow::Context;
use glob::glob;

const PID_FPR: f32 = 0.01;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("usage: pid <benchmark_dir>");
        return Ok(());
    };
    let bm_dir = Path::new(&args[1]);

    let benchmark = Benchmark::new(bm_dir.join("benchmark.tbl"))?;
    let data = PlotData::new(bm_dir.join("results"), &benchmark)?;

    Ok(())
}

struct Benchmark {
    target_cnt: usize,
    families: Vec<String>,
    targets_by_pid: HashMap<usize, Vec<String>>,
    intended_seq_queries_by_target: HashMap<String, String>,
}

impl Benchmark {
    fn new<P: AsRef<Path>>(tbl_path: P) -> anyhow::Result<Self> {
        let tbl_reader = BufReader::new(File::open(tbl_path)?);

        let mut target_cnt = 0;
        let mut families = HashSet::new();
        let mut targets_by_pid: HashMap<usize, Vec<String>> = HashMap::new();
        let mut queries_by_target = HashMap::new();

        for line in tbl_reader
            .lines()
            .map_while(Result::ok)
            .filter(|l| !l.starts_with('#'))
        {
            target_cnt += 1;
            let tokens: Vec<&str> = line.split_whitespace().collect();

            let pid = tokens[0]
                .strip_suffix('%')
                .context("pid entry missing % symbol")?
                .parse::<usize>()
                .context("")?;

            let family = tokens[1].to_string();
            let target = tokens[2].to_string();
            let query = tokens[3].to_string();

            families.insert(family);

            let targets = targets_by_pid.entry(pid).or_default();
            targets.push(target.clone());

            queries_by_target.insert(target, query);
        }

        Ok(Self {
            target_cnt,
            families: families.into_iter().collect::<Vec<_>>(),
            targets_by_pid,
            intended_seq_queries_by_target: queries_by_target,
        })
    }
}

pub struct PidData {
    name: String,
    points: Vec<(usize, f32)>,
}

impl PidData {
    fn new(tbl: &HitTable, bm: &Benchmark) -> Self {
        let decoy_cnt = (bm.target_cnt as f32 * PID_FPR).ceil() as usize;
        let mut positives_by_target: HashMap<String, Hit> = HashMap::new();
        let mut decoys_by_query: HashMap<String, Vec<Hit>> = HashMap::new();

        tbl.hits.iter().for_each(|h| {
            if h.target.starts_with("decoy") {
                let decoys = decoys_by_query.entry(h.query.clone()).or_default();
                decoys.push(h.clone());
            } else {
                let target_name = h.target.split('|').nth(1).unwrap().to_string();
                positives_by_target.insert(target_name, h.clone());
            }
        });

        let empty = vec![];
        let mut points = bm
            .targets_by_pid
            .iter()
            .map(|(&pid, targets)| {
                let total = targets.len() as f32;
                let found = targets
                    .iter()
                    .map(|target| match positives_by_target.get(target) {
                        Some(hit) => {
                            let decoys = decoys_by_query.get(&hit.query).unwrap_or(&empty);
                            // println!("{}", decoys.len());
                            // let e_value_thresh = decoys
                            //     .get(decoy_cnt + 1)
                            //     .expect("fewer decoys than FPR decoy count")
                            //     .e_value;

                            // if hit.e_value < e_value_thresh {
                            //     1.0
                            // } else {
                            //     0.0
                            // }
                            1.0
                        }
                        None => 0.0,
                    })
                    .sum::<f32>();

                (pid, found / total)
            })
            .collect::<Vec<_>>();

        points.sort_by(|a, b| a.0.cmp(&b.0));

        println!("target cnt: {}", bm.target_cnt);
        println!("decoy cnt: {}", decoy_cnt);
        println!("{}", tbl.name);
        points.iter().for_each(|p| println!("{p:?}"));

        Self {
            name: tbl.name.clone(),
            points,
        }
    }
}

pub struct RocData {
    name: String,
    pos: Vec<Hit>,
    decoys: Vec<Hit>,
}

impl RocData {
    fn new(tbl: &HitTable) -> Self {
        let (mut pos, mut decoys): (Vec<Hit>, Vec<Hit>) = tbl
            .hits
            .iter()
            .filter(|h| {
                let q = if h.query.contains("-consensus") {
                    h.query.split('-').next().unwrap()
                } else {
                    h.query.split('|').next().unwrap()
                };
                let t = h.target.split('|').next().unwrap();
                q == t || t.starts_with("decoy")
            })
            .cloned()
            .partition(|h| !h.target.starts_with("decoy"));

        pos.sort_by(|a, b| {
            a.e_value
                .partial_cmp(&b.e_value)
                .expect("NaN encountered in E-value sort")
        });

        decoys.sort_by(|a, b| {
            a.e_value
                .partial_cmp(&b.e_value)
                .expect("NaN encountered in E-value sort")
        });

        RocData {
            name: tbl.name.to_string(),
            pos,
            decoys,
        }
    }
}

pub struct PlotData {
    roc_data: Vec<RocData>,
    pid_data: Vec<PidData>,
    tools: Vec<String>,
}

impl PlotData {
    fn new<P: AsRef<Path>>(results_dir: P, bm: &Benchmark) -> anyhow::Result<Self> {
        let results_dir = PathBuf::from(results_dir.as_ref());
        let mut roc_data = vec![];
        let mut pid_data = vec![];
        let mut tools: HashSet<String> = HashSet::new();

        for path in glob(
            results_dir
                .join("*.tbl")
                .to_str()
                .context("invalid *.tbl glob")?,
        )?
        .filter_map(Result::ok)
        {
            let file = File::open(&path)?;
            let mut name = path
                .file_stem()
                .and_then(|s| s.to_str())
                .context("invalid path")?
                .splitn(2, '.');

            let tool = name.next().context("no tbl prefix")?;
            let search_type = name.next().context("no tbl search type")?;
            let name = format!("{tool} {search_type}");

            tools.insert(tool.to_string());

            let tbl = match tool {
                "hmmer" => HitTable::parse::<_, HmmerTable>(file, &name),
                "nail" => HitTable::parse::<_, NailTable>(file, &name),
                _ => HitTable::parse::<_, BlastTable>(file, &name),
            }?;

            match search_type {
                "cons" | "prf" => {
                    roc_data.push(RocData::new(&tbl));
                    pid_data.push(PidData::new(&tbl, bm));
                }
                "fam" => {
                    // what: filter the hit list such that for each (fam, target)
                    //       pair, retain only the best (lowest E-value) match
                    //
                    //  why: we want to smash family pairwise search results into
                    //       one set of hits per family so that we don't overcount
                    //       positives or decoys
                    //
                    //  how: build a hash that maps (fam, target) -> hit, where a hit
                    //       replaces an existing entry if it has a better E-value
                    let mut best_hit_by_fam_and_target: HashMap<(String, String), Hit> =
                        HashMap::new();

                    tbl.hits
                        .into_iter()
                        .map(|h| (h.query.split('|').next().unwrap().to_string(), h))
                        .for_each(|(q_fam, hit)| {
                            let key = (q_fam, hit.target.to_string());

                            match best_hit_by_fam_and_target.get(&key) {
                                Some(existing) => {
                                    if hit.e_value < existing.e_value {
                                        best_hit_by_fam_and_target.insert(key, hit);
                                    }
                                }
                                None => {
                                    best_hit_by_fam_and_target.insert(key, hit);
                                }
                            }
                        });

                    let tbl = HitTable {
                        name,
                        hits: best_hit_by_fam_and_target.into_values().collect(),
                    };

                    roc_data.push(RocData::new(&tbl));
                    pid_data.push(PidData::new(&tbl, bm));
                }
                "seq" => {
                    // what: filter out inter-family hits for non-intended
                    //       pairs of target/query sequences
                    //
                    //  why:
                    let intended_pair_hits: Vec<Hit> = tbl
                        .hits
                        .into_iter()
                        .filter(|h| {
                            if h.target.starts_with("decoy") {
                                return true;
                            }

                            let hit_target = h.target.split('|').nth(1).unwrap();
                            let hit_query = h.query.split('|').nth(1).unwrap();

                            let paired_query = bm
                                .intended_seq_queries_by_target
                                .get(hit_target)
                                .unwrap_or_else(|| panic!("{}", hit_target));

                            hit_query == *paired_query
                        })
                        .collect();

                    let tbl = HitTable {
                        name,
                        hits: intended_pair_hits,
                    };

                    // roc_data.push(RocData::new(&tbl));
                    pid_data.push(PidData::new(&tbl, bm));
                }
                _ => println!("unexpected search type flag: {search_type}"),
            };
        }

        Ok(Self {
            roc_data,
            pid_data,
            tools: tools.into_iter().collect(),
        })
    }
}
