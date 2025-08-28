use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    env,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use bioio::tbl::{parse_blast_tbl, parse_hmmer_domtbl, parse_nail_tbl, Hit};

use anyhow::anyhow;
use plotters::prelude::*;

const N_BINS: usize = 41;
const BIN_START: usize = 10;
const BIN_MAX: usize = N_BINS + BIN_START - 1;

const X_MIN: usize = 10;
const X_MAX: usize = 30;

const FPR: f32 = 0.001;

const TOL_PURPLE: &str = "#332288";
const TOL_GREEN: &str = "#117733";
const TOL_TEAL: &str = "#44AA99";
const TOL_BLUE: &str = "#88CCEE";
const TOL_YELLOW: &str = "#DDCC77";
const TOL_PINK: &str = "#CC6677";
const TOL_MAGENTA: &str = "#AA4499";
const TOL_OLIVE: &str = "#999933";

const COLORS: [&str; 8] = [
    TOL_PURPLE,
    TOL_GREEN,
    TOL_TEAL,
    TOL_BLUE,
    TOL_YELLOW,
    TOL_PINK,
    TOL_MAGENTA,
    TOL_OLIVE,
];

struct TargetInfo<'a> {
    family: &'a str,
    id: usize,
    bin: usize,
    start: usize,
    end: usize,
}

fn extract_target_info(hit: &Hit) -> anyhow::Result<TargetInfo> {
    //>Exo_endo_phos|83-289|10%:0 domain: O13348_MAGGR/17-223
    let mut tokens = hit.target.split('|');

    let family = tokens.next().ok_or(anyhow!("no family"))?;

    let mut range = tokens.next().ok_or(anyhow!("no range"))?.split('-');
    let start = range.next().ok_or(anyhow!("no start"))?.parse()?;
    let end = range.next().ok_or(anyhow!("no end"))?.parse()?;

    let mut binfo = tokens.next().ok_or(anyhow!("no bin"))?.split(':');
    let bin = binfo
        .next()
        .ok_or(anyhow!("no bin"))?
        .strip_suffix('%')
        .ok_or(anyhow!("no %"))?
        .parse()?;
    let id = binfo.next().ok_or(anyhow!("no id"))?.parse()?;

    Ok(TargetInfo {
        family,
        id,
        bin,
        start,
        end,
    })
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("usage: pid <benchmark_dir>");
        return Ok(());
    };

    let bm_dir = Path::new(&args[1]);
    let results_dir = bm_dir.join("results");

    let tbl = BufReader::new(File::open(bm_dir.join("benchmark.tbl"))?);

    let mut target_cnt_by_bin = [0usize; BIN_MAX + 1];
    let mut queries = HashSet::new();

    for line in tbl
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('#'))
    {
        let mut tokens = line.split_whitespace().next().unwrap().split('|');

        queries.insert(
            tokens
                .next()
                .expect("no query in benchmark table entry")
                .to_string(),
        );

        tokens
            .next_back()
            .expect("no bin in benchmark table entry")
            .split('%')
            .next()
            .expect("failed to split bin from id")
            .parse()
            .map(|b: usize| target_cnt_by_bin[b] += 1)?;
    }

    let queries = queries.iter().collect::<Vec<_>>();
    println!("queries found: {}", queries.len());

    println!("target bin distribution:");
    (BIN_START..=BIN_MAX)
        .collect::<Vec<usize>>()
        .chunks(5)
        .for_each(|bins| {
            bins.iter()
                .for_each(|b| print!("{b}%: {} | ", target_cnt_by_bin[*b]));
            println!();
        });

    let decoy_cnt = (FPR * queries.len() as f32) as usize;

    let mut hits = Vec::new();
    let name_fn = |p: PathBuf| Some((p.clone(), p.file_stem()?.to_str()?.to_string()));

    for (path, name) in glob::glob(results_dir.join("hmmer*.domtbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.push((name, parse_hmmer_domtbl(File::open(path)?)?));
    }

    for (path, name) in glob::glob(results_dir.join("mmseqs*.tbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.push((name, parse_blast_tbl(File::open(path)?)?));
    }

    for (path, name) in glob::glob(results_dir.join("nail*.tbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.push((name, parse_nail_tbl(File::open(path)?)?));
    }

    let root = SVGBackend::new("plot.svg", (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("CAPTION", ("sans-serif", 20))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(X_MIN..X_MAX, 0.0f32..1.0)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .x_label_formatter(&|x| format!("{}%", X_MAX - x + X_MIN))
        .draw()?;

    let filtered_hits = hits
        .into_iter()
        .map(|(name, hit_list)| {
            // what: filter the hit list such that for each (fam, target)
            //       pair, retain only the best (lowest E-value) match
            //
            //       also, filter out family mismatch hits while we're at it
            //
            //  why: since we're benchmarking pairwise-family search
            //
            let mut hit_hash: HashMap<(String, String), Hit> = HashMap::new();
            hit_list
                .into_iter()
                .map(|h| {
                    (
                        h.query.split('|').next().unwrap().to_string(),
                        h.target.split('|').next().unwrap().to_string(),
                        h,
                    )
                })
                // include only same-family hits & decoys
                .filter(|(q_fam, t_fam, _)| q_fam == t_fam || t_fam.starts_with("decoy"))
                .for_each(|(q_fam, _, hit)| {
                    let key = (q_fam, hit.target.to_string());

                    match hit_hash.entry(key) {
                        Entry::Occupied(mut entry) if hit.e_value < entry.get().e_value => {
                            entry.insert(hit);
                        }
                        Entry::Vacant(v) => {
                            v.insert(hit);
                        }
                        _ => {}
                    }
                });

            let mut filtered_hit_list: Vec<Hit> = hit_hash.into_values().collect();
            filtered_hit_list.sort_by(|a, b| {
                a.e_value
                    .partial_cmp(&b.e_value)
                    .expect("NaN E-value encountered")
            });

            (name, filtered_hit_list)
        })
        .collect::<Vec<_>>();

    let mut point_data = filtered_hits
        .into_iter()
        .map(|(name, hit_list)| {
            let mut num_hits_by_bin = vec![0usize; BIN_MAX + 1];

            let mut decoys_found = 0;
            for hit in hit_list.iter() {
                match hit.target.starts_with("decoy") {
                    true => {
                        decoys_found += 1;
                        if decoys_found >= decoy_cnt {
                            break;
                        }
                    }
                    false => {
                        let info = extract_target_info(hit).expect("failed to extract target info");
                        num_hits_by_bin[info.bin] += 1;
                    }
                }
            }

            let recall_by_bin: Vec<f32> = num_hits_by_bin
                .into_iter()
                .enumerate()
                .map(|(b, n)| {
                    let tot = target_cnt_by_bin[b];
                    if tot != 0 {
                        n as f32 / tot as f32
                    } else {
                        0.0
                    }
                })
                .collect();

            let points = recall_by_bin
                .into_iter()
                .enumerate()
                .skip(BIN_START)
                .map(|(b, r)| (X_MAX - b + X_MIN, r))
                .collect::<Vec<_>>();

            (name, points)
        })
        .collect::<Vec<_>>();

    point_data.sort_by(|a, b| {
        let s1: f32 = a.1.iter().map(|(_, y)| y).sum();
        let s2: f32 = b.1.iter().map(|(_, y)| y).sum();

        s1.partial_cmp(&s2).expect("NaN encountered")
    });

    let hex_to_rgb = |hex: &str| -> RGBColor {
        let hex = hex.trim_start_matches('#');
        let r = u8::from_str_radix(&hex[0..2], 16).unwrap();
        let g = u8::from_str_radix(&hex[2..4], 16).unwrap();
        let b = u8::from_str_radix(&hex[4..6], 16).unwrap();
        RGBColor(r, g, b)
    };

    for ((name, points), color) in point_data.iter().zip(COLORS) {
        let color = hex_to_rgb(color);
        chart
            .draw_series(LineSeries::new(points.clone(), color.stroke_width(3)))?
            .label(name)
            .legend(move |(x, y)| Rectangle::new([(x, y - 5), (x + 10, y + 5)], color.filled()));

        chart.draw_series(
            points
                .iter()
                .map(|(x, y)| Circle::new((*x, *y), 4, color.filled())),
        )?;

        chart.draw_series(
            points
                .iter()
                .map(|(x, y)| Circle::new((*x, *y), 1, WHITE.filled())),
        )?;
    }

    chart.configure_series_labels().border_style(BLACK).draw()?;

    Ok(())
}
