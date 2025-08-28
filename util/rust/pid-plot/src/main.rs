use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    env,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use bioio::tbl::{BlastTable, Hit, HitTable, HmmerDomTable, NailTable};

use anyhow::{anyhow, bail, Context};
use glob::glob;
use plotters::{element::DashedPathElement, prelude::*};

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

fn hex_to_rgb(hex: &str) -> RGBColor {
    let hex = hex.trim_start_matches('#');
    let r = u8::from_str_radix(&hex[0..2], 16).unwrap();
    let g = u8::from_str_radix(&hex[2..4], 16).unwrap();
    let b = u8::from_str_radix(&hex[4..6], 16).unwrap();
    RGBColor(r, g, b)
}

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

    let mut tables: Vec<HitTable> = vec![];
    let mut tools: HashSet<String> = HashSet::new();

    for path in glob(
        results_dir
            .join("*tbl")
            .to_str()
            .context("invalid *tbl glob")?,
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
        let suffix = name.next().context("no tbl suffix")?;
        let name = format!("{tool} {suffix}");

        tools.insert(tool.to_string());

        let tbl = match tool {
            "hmmer" => HitTable::parse::<_, HmmerDomTable>(file, &name),
            "nail" => HitTable::parse::<_, NailTable>(file, &name),
            _ => HitTable::parse::<_, BlastTable>(file, &name),
        }?;

        tables.push(tbl);
    }

    if tools.len() > COLORS.len() {
        bail!("not enough colors");
    }

    let mut tools: Vec<String> = tools.into_iter().collect();
    tools.sort();

    let color_by_tool: HashMap<String, RGBColor> = tools
        .into_iter()
        .zip(COLORS)
        .map(|(t, c)| (t, hex_to_rgb(c)))
        .collect();

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

    let filtered_tables: Vec<HitTable> = tables
        .into_iter()
        .map(|tbl| {
            // what: filter the hit list such that for each (fam, target)
            //       pair, retain only the best (lowest E-value) match
            //
            //       also, filter out family mismatch hits while we're at it
            //
            //  why: since we're benchmarking pairwise-family search
            //
            //  how: build a hash that maps (fam, target) -> hit, where a hit
            //       replaces an existing entry if it has a better E-value
            let mut hits_by_pair: HashMap<(String, String), Hit> = HashMap::new();
            tbl.hits
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

                    match hits_by_pair.get(&key) {
                        Some(existing) => {
                            if hit.e_value < existing.e_value {
                                hits_by_pair.insert(key, hit);
                            }
                        }
                        None => {
                            hits_by_pair.insert(key, hit);
                        }
                    }
                });

            let mut filtered_hits: Vec<Hit> = hits_by_pair.into_values().collect();
            filtered_hits.sort_by(|a, b| {
                a.e_value
                    .partial_cmp(&b.e_value)
                    .expect("NaN E-value encountered")
            });

            HitTable {
                name: tbl.name,
                hits: filtered_hits,
            }
        })
        .collect();

    let mut point_data = filtered_tables
        .into_iter()
        .map(|tbl| {
            let mut num_hits_by_bin = vec![0usize; BIN_MAX + 1];

            let mut decoys_found = 0;
            for hit in tbl.hits.iter() {
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

            (tbl.name, points)
        })
        .collect::<Vec<_>>();

    point_data.sort_by(|a, b| {
        let s1: f32 = a.1.iter().map(|(_, y)| y).sum();
        let s2: f32 = b.1.iter().map(|(_, y)| y).sum();

        s1.partial_cmp(&s2).expect("NaN encountered")
    });

    for (name, points) in point_data.iter() {
        let mut tokens = name.split_whitespace();
        let tool = tokens.next().context("no name")?;
        let search_type = tokens.next().context("no search type")?;
        let color = color_by_tool.get(tool).context("no color for tool")?;

        match search_type {
            "seq" => {
                chart
                    .draw_series(DashedLineSeries::new(
                        points.clone(),
                        5,
                        3,
                        color.stroke_width(3),
                    ))?
                    .label(name)
                    .legend(move |(x, y)| {
                        DashedPathElement::new(
                            vec![(x, y), (x + 21, y)],
                            5,
                            3,
                            color.stroke_width(3),
                        )
                    });
            }
            "prf" => {
                chart
                    .draw_series(LineSeries::new(points.clone(), color.stroke_width(3)))?
                    .label(name)
                    .legend(move |(x, y)| {
                        PathElement::new(vec![(x, y), (x + 21, y)], color.stroke_width(3))
                    });
            }
            _ => bail!("bad search type"),
        }

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
