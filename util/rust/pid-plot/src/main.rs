use std::{
    collections::{HashMap, HashSet},
    env,
    fs::File,
    io::{BufRead, BufReader},
    iter::once,
    path::{Path, PathBuf},
};

use bioio::tbl::{BlastTable, Hit, HitTable, HmmerTable, NailTable};

use anyhow::{anyhow, bail, Context};
use glob::glob;
use indexmap::IndexMap;
use once_cell::sync::Lazy;

use plotters::prelude::*;

const FPR: f32 = 0.01;

const PID_MAX: usize = 50;

const BIN_STEP: usize = 2;
const BIN_MIN: usize = 10;
const BIN_MAX: usize = 25;
const BIN_CNT: usize = (BIN_MAX - BIN_MIN) / BIN_STEP;

const X_MAX: usize = BIN_CNT + 2;

pub static PID_TO_X: Lazy<IndexMap<usize, usize>> = Lazy::new(|| {
    let mut m = IndexMap::new();
    (BIN_MIN..=BIN_MAX)
        .rev()
        .collect::<Vec<_>>()
        .chunks(BIN_STEP)
        .enumerate()
        .for_each(|(x, pids)| {
            for &pid in pids {
                m.insert(pid, x);
            }
        });
    m
});

pub static PID_RANGE_BY_X: Lazy<Vec<String>> = Lazy::new(|| {
    let mut v = vec![String::new(); X_MAX];
    (BIN_MIN..=BIN_MAX)
        .rev()
        .collect::<Vec<_>>()
        .chunks(BIN_STEP)
        .enumerate()
        .for_each(|(x, pids)| {
            v[x] = if pids.len() == 1 {
                format!("{}%", pids[0])
            } else {
                format!(
                    "{}-{}%",
                    pids.first().expect("empty pids in PID_RANGE_BY_X init"),
                    *pids.last().expect("empty pids in PID_RANGE_BY_X init")
                )
            };
        });
    v
});

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
    TOL_BLUE,
    TOL_YELLOW,
    TOL_PINK,
    TOL_TEAL,
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

fn extract_target_pid(name: &str) -> anyhow::Result<usize> {
    //>Exo_endo_phos|83-289|10%:0 domain: O13348_MAGGR/17-223
    Ok(name
        .split('|')
        .next_back()
        .ok_or(anyhow!("no bin"))?
        .split(':')
        .next()
        .ok_or(anyhow!("no bin"))?
        .strip_suffix('%')
        .ok_or(anyhow!("no %"))?
        .parse()?)
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("usage: pid <benchmark_dir>");
        return Ok(());
    };
    let bm_dir = Path::new(&args[1]);

    let benchmark = Benchmark::new(bm_dir.join("benchmark.tbl"))?;
    let tables = Tables::new(bm_dir.join("results"))?;

    lollipop(&benchmark, &tables)?;

    Ok(())
}

struct Benchmark {
    queries: Vec<String>,
    target_cnt_by_pid: [usize; PID_MAX + 1],
}

impl Benchmark {
    fn new<P: AsRef<Path>>(tbl_path: P) -> anyhow::Result<Self> {
        let tbl_reader = BufReader::new(File::open(tbl_path)?);

        let mut target_cnt_by_pid = [0usize; PID_MAX + 1];
        let mut queries = HashSet::new();

        for line in tbl_reader
            .lines()
            .map_while(Result::ok)
            .filter(|l| !l.starts_with('#'))
        {
            let mut tokens = line.split_whitespace();

            let pid = tokens
                .next()
                .context("no bin in benchmark table entry")?
                .split('%')
                .next()
                .context("failed to split bin from id")?
                .parse::<usize>()
                .context("failed to parse pid")?;

            target_cnt_by_pid[pid] += 1;

            queries.insert(
                tokens
                    .next_back()
                    .context("no query in benchmark table entry")?
                    .to_string(),
            );
        }

        Ok(Self {
            queries: queries.into_iter().collect::<Vec<_>>(),
            target_cnt_by_pid,
        })
    }
}

pub struct Tables {
    tables: Vec<HitTable>,
    tools: HashSet<String>,
}

impl Tables {
    fn new<P: AsRef<Path>>(results_dir: P) -> anyhow::Result<Self> {
        let results_dir = PathBuf::from(results_dir.as_ref());
        let mut tables: Vec<HitTable> = vec![];
        let mut tools: HashSet<String> = HashSet::new();

        for path in glob(
            results_dir
                .join("*.tbl")
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
                "hmmer" => HitTable::parse::<_, HmmerTable>(file, &name),
                "nail" => HitTable::parse::<_, NailTable>(file, &name),
                _ => HitTable::parse::<_, BlastTable>(file, &name),
            }?;

            tables.push(tbl);
        }

        let tables: Vec<HitTable> = tables
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
                        if h.query.contains("consensus") {
                            (
                                h.query.split('-').next().unwrap().to_string(),
                                h.target.split('|').next().unwrap().to_string(),
                                h,
                            )
                        } else {
                            (
                                h.query.split('|').next().unwrap().to_string(),
                                h.target.split('|').next().unwrap().to_string(),
                                h,
                            )
                        }
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

        Ok(Self { tables, tools })
    }
}

fn lollipop(benchmark: &Benchmark, tables: &Tables) -> anyhow::Result<()> {
    let decoy_cnt = (FPR * benchmark.queries.len() as f32).ceil() as usize;

    let mut target_cnt_by_x = [0usize; X_MAX];
    benchmark
        .target_cnt_by_pid
        .iter()
        .enumerate()
        .for_each(|(pid, cnt)| {
            if (BIN_MIN..=BIN_MAX).contains(&pid) {
                let x = PID_TO_X.get(&pid).unwrap();
                target_cnt_by_x[*x] += cnt;
            }
        });

    let mut point_data = tables
        .tables
        .iter()
        .map(|tbl| {
            let mut num_hits_by_x = vec![0usize; X_MAX - 1];

            let mut decoys_found = 0;
            for hit in tbl.hits.iter() {
                match hit.target.starts_with("decoy") {
                    true => {
                        decoys_found += 1;
                        if decoys_found >= decoy_cnt {
                            println!("{} {:.2} {:.2e}", tbl.name, hit.score, hit.e_value);
                            break;
                        }
                    }
                    false => {
                        let pid =
                            extract_target_pid(&hit.target).expect("failed to extract target info");

                        if (BIN_MIN..=BIN_MAX).contains(&pid) {
                            let x = PID_TO_X
                                .get(&pid)
                                .unwrap_or_else(|| panic!("no x for pid: {pid}",));
                            num_hits_by_x[*x] += 1;
                        }
                    }
                }
            }

            let recall_by_x: Vec<f32> = num_hits_by_x
                .into_iter()
                .enumerate()
                .map(|(x, n)| {
                    let tot = target_cnt_by_x[x];
                    if tot != 0 {
                        n as f32 / tot as f32
                    } else {
                        0.0
                    }
                })
                .collect();

            let points = recall_by_x
                .into_iter()
                .enumerate()
                .map(|(x, recall)| (x as f32, recall))
                .collect::<Vec<_>>();

            (tbl.name.clone(), points)
        })
        .collect::<Vec<_>>();

    point_data.sort_by(|a, b| {
        let s1: f32 = a.1.iter().map(|(_, y)| y).sum();
        let s2: f32 = b.1.iter().map(|(_, y)| y).sum();

        s1.partial_cmp(&s2).expect("NaN encountered")
    });

    point_data.reverse();
    if tables.tools.len() > COLORS.len() {
        bail!("not enough colors");
    }

    let mut tools: Vec<String> = tables.tools.iter().cloned().collect();
    tools.sort();

    let color_by_tool: HashMap<String, RGBColor> = tools
        .into_iter()
        .zip(COLORS)
        .map(|(t, c)| (t, hex_to_rgb(c)))
        .collect();

    let width = 600;
    let height = 500;
    let root = SVGBackend::new("pid.svg", (width, height)).into_drawing_area();
    root.fill(&WHITE)?;
    let (top, bottom) = root.split_vertically((height / 4) * 3);

    let x_label_fmt = |x: &f32| -> String { PID_RANGE_BY_X.get(*x as usize).unwrap().clone() };

    let black = RGBColor(12, 4, 4);
    let grey = &full_palette::GREY_800;
    let x_range = 0f32..(X_MAX - 1) as f32;
    let axis_style = grey.stroke_width(2);
    let axis_desc_style = ("sans-serif", 12, &black);
    let label_style = ("Arial", 8, grey);

    let mut recall_chart = ChartBuilder::on(&top)
        .caption(
            "Recall by decreasing % pairwise identity",
            ("sans-serif", 20),
        )
        .margin(20)
        .margin_bottom(0)
        .x_label_area_size(18)
        .top_x_label_area_size(40)
        .y_label_area_size(40)
        .right_y_label_area_size(40)
        .build_cartesian_2d(x_range.clone(), 0.0f32..1.0)?
        .set_secondary_coord(x_range.clone(), 0.0f32..1.0);

    recall_chart
        .configure_mesh()
        .disable_mesh()
        .x_labels(X_MAX)
        .x_label_formatter(&x_label_fmt)
        .y_desc(format!("fraction of TP found at {FPR} FP per query"))
        .axis_style(axis_style)
        .axis_desc_style(axis_desc_style)
        .label_style(label_style)
        .draw()?;

    recall_chart
        .configure_secondary_axes()
        .x_labels(X_MAX)
        .x_label_formatter(&x_label_fmt)
        .y_desc("")
        .axis_style(axis_style)
        .axis_desc_style(axis_desc_style)
        .label_style(label_style)
        .draw()?;

    recall_chart.draw_series(DashedLineSeries::new(
        [(0.0, 0.5), (X_MAX as f32, 0.5)],
        10,
        7,
        BLACK.mix(0.50).stroke_width(1),
    ))?;

    let max_cnt = (*target_cnt_by_x.iter().max().unwrap() as f32 / 100.0).round() * 100.0;

    let mut bin_chart = ChartBuilder::on(&bottom)
        .margin(20)
        .margin_top(0)
        .x_label_area_size(10)
        .top_x_label_area_size(10)
        .y_label_area_size(40)
        .right_y_label_area_size(40)
        .build_cartesian_2d(x_range.clone(), 0.0f32..max_cnt)?
        .set_secondary_coord(x_range.clone(), 0.0f32..max_cnt);

    bin_chart
        .configure_mesh()
        .disable_mesh()
        .x_labels(X_MAX)
        .x_label_formatter(&x_label_fmt)
        .y_labels(6)
        .y_label_formatter(&|y| format!("{}", max_cnt - y))
        .y_desc("# of sequence pairs")
        .axis_style(axis_style)
        .axis_desc_style(axis_desc_style)
        .label_style(label_style)
        .draw()?;

    bin_chart
        .configure_secondary_axes()
        .x_labels(X_MAX)
        .x_label_formatter(&|_| "".to_string())
        .y_labels(6)
        .y_label_formatter(&|y| format!("{}", max_cnt - y))
        .y_desc("")
        .axis_style(axis_style)
        .axis_desc_style(axis_desc_style)
        .label_style(label_style)
        .draw()?;

    let bin_cnt_points = target_cnt_by_x
        .iter()
        .enumerate()
        .map(|(x, cnt)| (x as f32, max_cnt - *cnt as f32))
        .collect::<Vec<_>>();

    let bin_chart_color = hex_to_rgb(TOL_MAGENTA);

    bin_chart.draw_series([Polygon::new(
        once((0.0f32, max_cnt))
            .chain(bin_cnt_points.iter().cloned())
            .chain(once((PID_MAX as f32, max_cnt)))
            .collect::<Vec<_>>(),
        bin_chart_color.mix(0.60),
    )])?;

    bin_chart.draw_series(LineSeries::new(
        bin_cnt_points.clone(),
        bin_chart_color.stroke_width(2),
    ))?;

    bin_chart.draw_series(
        bin_cnt_points
            .iter()
            .map(|(x, y)| Circle::new((*x, *y), 3, bin_chart_color.filled())),
    )?;

    bin_chart.draw_series(
        bin_cnt_points
            .iter()
            .map(|(x, y)| Circle::new((*x, *y), 1, WHITE.filled())),
    )?;

    for (i, (name, points)) in point_data.iter().enumerate() {
        let mut tokens = name.split_whitespace();
        let tool = tokens.next().context("no name")?;
        let search_type = tokens.next().context("no search type")?;
        let color = color_by_tool.get(tool).context("no color for tool")?;

        let o1 = i as f32 / 15.0;
        recall_chart.draw_series(points.iter().map(|(x, y)| {
            PathElement::new(vec![(*x + o1, 0.0), (*x + o1, *y)], color.stroke_width(3))
        }))?;

        match search_type {
            "prf" => {
                recall_chart
                    .draw_series(
                        points
                            .iter()
                            .map(|(x, y)| Circle::new((*x + o1, *y), 2, color.filled())),
                    )?
                    .label(name)
                    .legend(move |(x, y)| {
                        EmptyElement::at((0, 0))
                            + PathElement::new(vec![(x, y), (x + 21, y)], color.stroke_width(3))
                            + Circle::new((x + 21, y), 2, color.filled())
                    });
            }
            "pair" => {
                recall_chart
                    .draw_series(
                        points
                            .iter()
                            .map(|(x, y)| Circle::new((*x + o1, *y), 2, color.filled())),
                    )?
                    .label(name)
                    .legend(move |(x, y)| {
                        EmptyElement::at((0, 0))
                            + PathElement::new(vec![(x, y), (x + 21, y)], color.stroke_width(3))
                            + Circle::new((x + 21, y), 2, color.filled())
                            + Circle::new((x + 21, y), 1, BLACK.filled())
                    });

                recall_chart.draw_series(
                    points
                        .iter()
                        .map(|(x, y)| Circle::new((*x + o1, *y), 1, BLACK.filled())),
                )?;
            }
            "cons" => {
                recall_chart
                    .draw_series(
                        points
                            .iter()
                            .map(|(x, y)| Circle::new((*x + o1, *y), 2, color.filled())),
                    )?
                    .label(name)
                    .legend(move |(x, y)| {
                        EmptyElement::at((0, 0))
                            + PathElement::new(vec![(x, y), (x + 21, y)], color.stroke_width(3))
                            + Circle::new((x + 21, y), 2, color.filled())
                            + Circle::new((x + 21, y), 1, WHITE.filled())
                    });

                recall_chart.draw_series(
                    points
                        .iter()
                        .map(|(x, y)| Circle::new((*x + o1, *y), 1, WHITE.filled())),
                )?;
            }
            _ => bail!("bad search type"),
        }
    }

    recall_chart
        .configure_series_labels()
        .border_style(BLACK)
        .margin(5)
        .background_style(WHITE.filled())
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    Ok(())
}
