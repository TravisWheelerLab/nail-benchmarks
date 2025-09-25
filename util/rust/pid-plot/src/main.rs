use std::{
    collections::{HashMap, HashSet},
    env,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, Write},
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

    data.write_pid(&mut std::io::stdout())?;
    Ok(())
}

struct BenchmarkEntry {
    pid: usize,
    target: String,
    query: String,
    family: String,
}

struct Benchmark {
    target_cnt: usize,
    entries: Vec<BenchmarkEntry>,
    idx_by_pid: HashMap<usize, Vec<usize>>,
    idx_by_target: HashMap<String, usize>,
    idx_by_query: HashMap<String, Vec<usize>>,
    idx_by_family: HashMap<String, Vec<usize>>,
}

impl Benchmark {
    fn new<P: AsRef<Path>>(tbl_path: P) -> anyhow::Result<Self> {
        let tbl_reader = BufReader::new(File::open(tbl_path)?);

        let mut target_cnt = 0;
        let mut entries = vec![];
        let mut idx_by_pid: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut idx_by_target: HashMap<String, usize> = HashMap::new();
        let mut idx_by_query: HashMap<String, Vec<usize>> = HashMap::new();
        let mut idx_by_family: HashMap<String, Vec<usize>> = HashMap::new();

        for line in tbl_reader
            .lines()
            .map_while(Result::ok)
            .filter(|l| !l.starts_with('#'))
        {
            let tokens: Vec<&str> = line.split_whitespace().collect();

            let pid = tokens[0]
                .strip_suffix('%')
                .context("pid entry missing % symbol")?
                .parse::<usize>()
                .context("")?;
            let family = tokens[1].to_string();
            let target = tokens[2].to_string();
            let query = tokens[3].to_string();

            let entry = BenchmarkEntry {
                pid,
                target: target.clone(),
                query: query.clone(),
                family: family.clone(),
            };

            idx_by_pid.entry(pid).or_default().push(target_cnt);
            idx_by_target.insert(target, target_cnt);
            idx_by_query.entry(query).or_default().push(target_cnt);
            idx_by_family.entry(family).or_default().push(target_cnt);

            entries.push(entry);
            target_cnt += 1;
        }

        Ok(Self {
            target_cnt,
            entries,
            idx_by_pid,
            idx_by_target,
            idx_by_query,
            idx_by_family,
        })
    }
}

impl Benchmark {
    fn entries_by_pid(&self, pid: usize) -> Vec<&BenchmarkEntry> {
        self.idx_by_pid
            .get(&pid)
            .unwrap()
            .iter()
            .map(|i| &self.entries[*i])
            .collect()
    }

    fn entry_by_target(&self, target: &str) -> &BenchmarkEntry {
        &self.entries[*self.idx_by_target.get(target).unwrap()]
    }

    fn entries_by_query(&self, query: &str) -> Vec<&BenchmarkEntry> {
        self.idx_by_query
            .get(query)
            .unwrap()
            .iter()
            .map(|i| &self.entries[*i])
            .collect()
    }

    fn entries_by_family(&self, family: &str) -> Vec<&BenchmarkEntry> {
        self.idx_by_family
            .get(family)
            .unwrap()
            .iter()
            .map(|i| &self.entries[*i])
            .collect()
    }
}

#[derive(PartialEq)]
enum SearchType {
    Profile,
    Consensus,
    Sequence,
}

impl Display for SearchType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            SearchType::Profile => "prf",
            SearchType::Consensus => "cons",
            SearchType::Sequence => "seq",
        };
        write!(f, "{s}")
    }
}

struct PidData {
    name: String,
    points: Vec<(usize, f32)>,
}

impl PidData {
    fn new(tbl: &HitTable, bm: &Benchmark, search_type: SearchType) -> Self {
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

        let mut complete_decoys: Vec<Hit> = bm
            .entries
            .iter()
            .flat_map(|entry| {
                let name = match search_type {
                    SearchType::Profile => &entry.family,
                    SearchType::Consensus => &format!("{}-consensus", entry.family),
                    SearchType::Sequence => &format!("{}|{}", entry.family, entry.query),
                };
                match decoys_by_query.get(name) {
                    Some(hits) => hits.clone(),
                    None => vec![],
                }
            })
            .collect();

        complete_decoys.sort_by(|a, b| {
            a.e_value
                .partial_cmp(&b.e_value)
                .expect("NaN encountered while sorting complete decoy list")
        });

        let n_decoys_allowed = (bm.target_cnt as f32 * PID_FPR).ceil() as usize;

        let e_value_thresh = complete_decoys
            .get(n_decoys_allowed + 1)
            .expect("not enough decoys to produce E-value threshold")
            .e_value;

        let mut points = bm
            .idx_by_pid
            .keys()
            .map(|pid| {
                (
                    pid,
                    bm.entries_by_pid(*pid)
                        .iter()
                        .map(|e| &e.target)
                        .collect::<Vec<_>>(),
                )
            })
            .map(|(&pid, targets)| {
                let total = targets.len() as f32;
                let found = targets
                    .iter()
                    .map(|target| match positives_by_target.get(*target) {
                        Some(hit) => {
                            if hit.e_value < e_value_thresh {
                                1.0
                            } else {
                                0.0
                            }
                        }
                        None => 0.0,
                    })
                    .sum::<f32>();

                (pid, found / total)
            })
            .collect::<Vec<_>>();

        points.sort_by(|a, b| a.0.cmp(&b.0));

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
            let search_type = match name.next().context("no tbl search type")? {
                "cons" => SearchType::Consensus,
                "prf" => SearchType::Profile,
                "seq" => SearchType::Sequence,
                _ => panic!("unknown search type"),
            };

            let name = format!("{tool} {search_type}");

            tools.insert(tool.to_string());

            let tbl = match tool {
                "hmmer" => HitTable::parse::<_, HmmerTable>(file, &name),
                "nail" => HitTable::parse::<_, NailTable>(file, &name),
                _ => HitTable::parse::<_, BlastTable>(file, &name),
            }?;

            match search_type {
                SearchType::Consensus | SearchType::Profile => {
                    roc_data.push(RocData::new(&tbl));
                    pid_data.push(PidData::new(&tbl, bm, search_type));
                }
                SearchType::Sequence => {
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
                            let paired_query = &bm.entry_by_target(hit_target).query;

                            hit_query == paired_query
                        })
                        .collect();

                    let tbl = HitTable {
                        name,
                        hits: intended_pair_hits,
                    };

                    pid_data.push(PidData::new(&tbl, bm, search_type));
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

    fn write_pid<W: Write>(&self, out: &mut W) -> anyhow::Result<()> {
        self.pid_data.iter().try_for_each(|d| {
            write!(out, "{}", d.name)?;
            d.points
                .iter()
                .try_for_each(|p| write!(out, ",({}, {:.3})", p.0, p.1))?;
            writeln!(out)
        })?;
        Ok(())
    }
}
