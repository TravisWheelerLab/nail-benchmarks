use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use bioio::tbl::{
    hmmer::{HmmerDomainTable, HmmerHit},
    BlastTable, Hit, HitTable, HmmerTable, NailTable,
};

pub fn p_value(score: f64, lambda: f64, tau: f64) -> f64 {
    (-lambda * (score - tau)).exp()
}

struct HmmGumbel {
    ga_full: f64,
    ga_dom: f64,
    tau: f64,
    lambda: f64,
}

impl HmmGumbel {
    pub fn p_value(&self, score: f64) -> f64 {
        (-self.lambda * (score - self.tau)).exp()
    }
}

#[derive(Debug)]
enum Score {
    None,
    Pass(f32),
    Filtered(f32),
}

#[derive(Debug)]
struct NailStats {
    query: String,
    target: String,
    cloud_score: Score,
    forward_score: Score,
}

fn nail_stats(
    stats_path: impl AsRef<Path>,
) -> anyhow::Result<HashMap<(String, String), NailStats>> {
    let reader = BufReader::new(File::open(stats_path.as_ref())?);

    let mut stats = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = &line[1..line.len()];
        let tokens = line
            .split(") (")
            .map(|t| t.split_whitespace().collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let query = tokens[0][0].to_string();
        let target = tokens[0][1].to_string();

        let sc = tokens[1][1]
            .strip_suffix("b")
            .unwrap()
            .parse::<f32>()
            .unwrap();

        let cloud_score = if tokens[1][0] == "P" {
            Score::Pass(sc)
        } else {
            Score::Filtered(sc)
        };

        let forward_score = if tokens[1][0] == "P" {
            let s = tokens[2][1]
                .strip_suffix("b")
                .unwrap()
                .parse::<f32>()
                .unwrap();
            if tokens[2][0] == "P" {
                Score::Pass(s)
            } else {
                Score::Filtered(s)
            }
        } else {
            Score::None
        };

        let s = NailStats {
            query,
            target,
            cloud_score,
            forward_score,
        };

        stats.insert((s.query.clone(), s.target.clone()), s);
    }

    Ok(stats)
}

fn ga_pvals(hmm_path: impl AsRef<Path>) -> anyhow::Result<HashMap<String, HmmGumbel>> {
    let reader = BufReader::new(File::open(hmm_path.as_ref())?);

    let mut names = vec![];
    let mut gathering_thresholds = vec![];
    let mut gumbels = vec![];

    for line in reader.lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix("NAME") {
            names.push(rest.split_whitespace().collect::<String>());
        }

        if let Some(rest) = line.strip_prefix("GA") {
            let x = rest
                .split_whitespace()
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>();

            gathering_thresholds.push((x[0], x[1]));
        }

        if let Some(rest) = line.strip_prefix("STATS LOCAL FORWARD") {
            let x = rest
                .split_whitespace()
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>();

            gumbels.push((x[0], x[1]))
        }
    }

    assert_eq!(names.len(), gathering_thresholds.len());
    assert_eq!(names.len(), gumbels.len());

    Ok(names
        .into_iter()
        .enumerate()
        .map(|(i, n)| {
            (
                n,
                HmmGumbel {
                    ga_full: gathering_thresholds[1].0,
                    ga_dom: gathering_thresholds[1].1,
                    tau: gumbels[i].0,
                    lambda: gumbels[i].1,
                },
            )
        })
        .collect())
}

#[derive(Default)]
struct Record {
    query: String,
    target: String,
    ga_score: f32,
    ga_p_value: f64,
    nail_cloud_score: Option<f32>,
    nail_cloud_p_value: Option<f64>,
    nail_score: Option<f32>,
    nail_p_value: Option<f64>,
    mmseqs_score: Option<f32>,
    mmseqs_p_value: Option<f64>,
    hmmer_score: Option<f32>,
    hmmer_p_value: Option<f64>,
    dom_score_sum: Option<f32>,
    dom_score_max: Option<f32>,
    dom_sig_cnt: Option<usize>,
    dom_scores: Option<Vec<f32>>,
}

impl Record {
    pub fn new(
        hmmer_hit: Option<&Hit>,
        hmmer_dom: Option<&HmmerHit>,
        nail_hit: Option<&Hit>,
        nail_stats: Option<&NailStats>,
        mmseqs_hit: Option<&Hit>,
        hmm: &HmmGumbel,
        z: f64,
    ) -> Self {
        let hit = hmmer_hit.or(nail_hit).or(mmseqs_hit).unwrap();
        let (query, target) = (hit.query.to_string(), hit.target.to_string());

        let ga_score = hmm.ga_full as f32;
        let ga_p_value = hmm.p_value(hmm.ga_full);

        let (nail_score, nail_p_value) = nail_hit
            .map(|h| (Some(h.score as f32), Some(h.e_value / z)))
            .unwrap_or_default();

        let (nail_cloud_score, nail_cloud_p_value) = nail_stats
            .map(|s| match s.cloud_score {
                Score::None => (None, None),
                Score::Pass(sc) | Score::Filtered(sc) => (Some(sc), Some(hmm.p_value(sc as f64))),
            })
            .unwrap_or_default();

        let (mmseqs_score, mmseqs_p_value) = mmseqs_hit
            .map(|h| (Some(h.score as f32), Some(h.e_value / z)))
            .unwrap_or_default();

        let (hmmer_score, hmmer_p_value) = hmmer_hit
            .map(|h| (Some(h.score as f32), Some(h.e_value / z)))
            .unwrap_or_default();

        let (dom_score_sum, dom_score_max, dom_sig_cnt, dom_scores) = hmmer_dom
            .map(|dom| {
                let mut dom_scores = dom.domains.iter().map(|d| d.score).collect::<Vec<_>>();

                let dom_score_sum = dom_scores.iter().sum();
                dom_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let dom_score_max = *dom_scores.last().unwrap();

                dom_scores.reverse();

                let dom_pct = dom_scores
                    .iter()
                    .map(|s| s / dom_score_max)
                    .collect::<Vec<_>>();
                let dom_sig_cnt = dom_pct.iter().filter(|p| **p >= 0.1).count();
                (
                    Some(dom_score_sum),
                    Some(dom_score_max),
                    Some(dom_sig_cnt),
                    Some(dom_scores),
                )
            })
            .unwrap_or_default();

        Self {
            query,
            target,
            ga_score,
            ga_p_value,
            nail_cloud_score,
            nail_cloud_p_value,
            nail_score,
            nail_p_value,
            mmseqs_score,
            mmseqs_p_value,
            hmmer_score,
            hmmer_p_value,
            dom_score_sum,
            dom_score_max,
            dom_sig_cnt,
            dom_scores,
        }
    }
}

fn score_fmt(f: Option<f32>, w: usize) -> String {
    match f {
        Some(s) => format!("{s:<W$.1}", W = w),
        None => format!("{:<W$}", "-", W = w),
    }
}

fn p_value_fmt(f: Option<f64>, w: usize) -> String {
    match f {
        Some(s) => format!("{s:<W$.2e}", W = w),
        None => format!("{:<W$}", "-", W = w),
    }
}

fn int_fmt(f: Option<usize>, w: usize) -> String {
    match f {
        Some(s) => format!("{s:<W$}", W = w),
        None => format!("{:<W$}", "-", W = w),
    }
}

const HEADER: [&str; 3] = [
    "| GA |    GA     |  nail  |    nail   |  nail  |   nail    | mmseqs |  mmseqs   | hmmer |  hmmer    | dom | dom | sig |",
    "| sc |  P-value  | cld sc |cld P-value|   sc   |  P-value  |   sc   |  P-value  |  sc   |  P-value  | sum | max | dom | dom scores",
    " ---- ----------- -------- ----------- -------- ----------- -------- ----------- ------- ----------- ----- ----- ----- ------------",
    ];

fn write_header(out: &mut impl Write, recs: &[Record]) -> anyhow::Result<()> {
    let q_max = recs.iter().map(|r| r.query.len()).max().unwrap() + 1;
    let t_max = recs.iter().map(|r| r.target.len()).max().unwrap() + 1;

    writeln!(
        out,
        "{:W1$} {:W2$} {}",
        " ",
        " ",
        HEADER[0],
        W1 = q_max,
        W2 = t_max
    )?;
    writeln!(
        out,
        "{:W1$} {:W2$} {}",
        "query",
        "target",
        HEADER[1],
        W1 = q_max,
        W2 = t_max,
    )?;
    writeln!(
        out,
        "{:W1$} {:W2$} {}",
        "-".repeat(q_max),
        "-".repeat(t_max),
        HEADER[2],
        W1 = q_max,
        W2 = t_max,
    )?;

    Ok(())
}

fn write_records(out: &mut impl Write, recs: &[Record], w: &[usize]) -> anyhow::Result<()> {
    let q_max = recs.iter().map(|r| r.query.len()).max().unwrap() + 1;
    let t_max = recs.iter().map(|r| r.target.len()).max().unwrap() + 1;

    for r in recs.iter() {
        writeln!(
            out,
            "{:W1$} {:W2$}  {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            r.query,
            r.target,
            score_fmt(Some(r.ga_score), w[0]),
            p_value_fmt(Some(r.ga_p_value), w[1]),
            score_fmt(r.nail_cloud_score, w[2]),
            p_value_fmt(r.nail_cloud_p_value, w[3]),
            score_fmt(r.nail_score, w[4]),
            p_value_fmt(r.nail_p_value, w[5]),
            score_fmt(r.mmseqs_score, w[6]),
            p_value_fmt(r.mmseqs_p_value, w[7]),
            score_fmt(r.hmmer_score, w[8]),
            p_value_fmt(r.hmmer_p_value, w[9]),
            score_fmt(r.dom_score_sum, w[10]),
            score_fmt(r.dom_score_max, w[11]),
            int_fmt(r.dom_sig_cnt, w[12]),
            match r.dom_scores {
                Some(ref v) => v
                    .iter()
                    .map(|f| format!("{f:.1}"))
                    .collect::<Vec<_>>()
                    .join(","),
                None => "-".to_string(),
            },
            W1 = q_max,
            W2 = t_max
        )?;
    }
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let z = 2_455_939.0;

    let hmm_path =
        PathBuf::from("/Users/jack/projects/nail-benchmarks/mgnify-pfam/benchmark-sm/query.hmm");

    let dir = PathBuf::from("/Users/jack/projects/nail-benchmarks/mgnify-pfam/1/");
    let hmmer_tbl = HitTable::from_path::<_, HmmerTable>(dir.join("hmmer.1.prf.tbl"))?;
    let hmmer_domtbl = HmmerDomainTable::from_path(dir.join("hmmer.1.prf.domtbl"))?;
    let nail_tbl = HitTable::from_path::<_, NailTable>(dir.join("nail.1.prf.tbl"))?;
    let nail_stats = nail_stats(dir.join("nail.1.prf.stats"))?;
    let mmseqs_tbl = HitTable::from_path::<_, BlastTable>(dir.join("mmseqs.1.prf.tbl"))?;

    let hmm_gumbels = ga_pvals(hmm_path)?;
    let mut queries = hmm_gumbels.keys().cloned().collect::<Vec<_>>();
    queries.sort();

    let mut hits_by_query = HashMap::new();

    fn map_fn(
        tbl: HitTable,
        hits: &mut HashMap<String, Vec<String>>,
    ) -> HashMap<(String, String), Hit> {
        tbl.hits
            .into_iter()
            .map(|h| {
                let vec = hits.entry(h.query.clone()).or_default();

                if !vec.contains(&h.target) {
                    vec.push(h.target.clone());
                }

                ((h.query.clone(), h.target.clone()), h)
            })
            .collect::<HashMap<_, _>>()
    }

    let hmmer_map = map_fn(hmmer_tbl, &mut hits_by_query);
    let nail_map = map_fn(nail_tbl, &mut hits_by_query);
    let mmseqs_map = map_fn(mmseqs_tbl, &mut hits_by_query);

    let hits: Vec<(String, String)> = queries
        .into_iter()
        .flat_map(|q| {
            hits_by_query
                .remove(&q)
                .unwrap_or_default()
                .into_iter()
                .map(|t| (q.clone(), t.clone()))
                .collect::<Vec<(String, String)>>()
        })
        .collect();

    let mut header_written = false;
    let widths = HEADER[2]
        .split_whitespace()
        .map(|s| s.len())
        .collect::<Vec<usize>>();
    let mut recs: Vec<Record> = vec![];
    let mut out = BufWriter::new(File::create("out.tbl")?);
    for k in hits {
        let hmm = hmm_gumbels.get(&k.0).unwrap();

        let hmmer_hit = hmmer_map.get(&k);
        let nail_hit = nail_map.get(&k);
        let mmseqs_hit = mmseqs_map.get(&k);

        let hmmer_passed = hmmer_hit.map(|h| h.score >= hmm.ga_full).unwrap_or(false);
        let nail_passed = nail_hit.map(|h| h.score >= hmm.ga_full).unwrap_or(false);
        let mmseqs_passed = mmseqs_hit
            .map(|h| h.e_value / z <= hmm.p_value(hmm.ga_full))
            .unwrap_or(false);

        if !(hmmer_passed || nail_passed || mmseqs_passed) {
            continue;
        }

        let hmmer_dom = hmmer_domtbl.hits.get(&k);
        let nail_stats = nail_stats.get(&k);

        let rec = Record::new(
            hmmer_hit, hmmer_dom, nail_hit, nail_stats, mmseqs_hit, hmm, z,
        );

        if !recs.is_empty() && rec.query != *recs.last().unwrap().query {
            if !header_written {
                write_header(&mut out, &recs)?;
                header_written = true;
            }

            write_records(&mut out, &recs, &widths)?;
            recs.clear();
        }

        recs.push(rec);
    }
    write_records(&mut out, &recs, &widths)?;

    Ok(())
}
