use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    env,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use bioio::tools::{hmmer, mmseqs, nail, Hit};

use anyhow::{anyhow, bail};

use plotly::{common::Mode, ImageFormat, Layout, Plot, Scatter};

const N_BINS: usize = 41;
const BIN_START: usize = 10;
const MAX_BIN: usize = N_BINS + BIN_START - 1;
const FPR: f32 = 0.01;

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

fn template() -> plotly::layout::Template {
    use plotly::{
        common::{ColorBar, Font, Title},
        layout::{Axis, AxisRange, ColorAxis, LayoutTemplate, Template, TicksDirection},
    };
    let layout_template = LayoutTemplate::new()
        .color_axis(ColorAxis::new().color_bar(ColorBar::new().outline_width(0)))
        .colorway(COLORS.to_vec())
        .font(Font::new().color("#2a3f5f"))
        .paper_background_color("#FFFFFF")
        .plot_background_color("#FFFFFF")
        .title(Title::new().x(0.05))
        .x_axis(
            Axis::new()
                .auto_margin(true)
                .show_line(true)
                .line_color("black")
                .line_width(2)
                .ticks(TicksDirection::Outside)
                .tick_width(2)
                .tick_color("black")
                .mirror(true)
                .grid_color("#e5e5e5")
                .zero_line_color("#e5e5e5")
                .zero_line_width(1)
                .tick_suffix("%")
                .range(AxisRange::new(50, 0)),
        )
        .y_axis(
            Axis::new()
                .auto_margin(true)
                .show_line(true)
                .line_color("black")
                .line_width(2)
                .ticks(TicksDirection::Outside)
                .tick_width(2)
                .tick_color("black")
                .mirror(true)
                .grid_color("#e5e5e5")
                .zero_line_color("#e5e5e5")
                .zero_line_width(1)
                .range(AxisRange::new(0, 1.0)),
        );
    Template::new().layout(layout_template)
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
        println!("usage: pid <benchmark_dir> [figures]");
        return Ok(());
    };

    let bm_dir = Path::new(&args[1]);
    let results_dir = bm_dir.join("results");

    let fig_dir = if args.len() > 2 {
        Path::new(&args[2])
    } else {
        Path::new("./figures/")
    };

    let tbl = BufReader::new(File::open(bm_dir.join("benchmark.tbl"))?);

    let mut target_cnt_by_bin = [0usize; MAX_BIN + 1];
    let mut queries = HashSet::new();

    for line in tbl
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('#'))
    {
        let tokens: Vec<&str> = line.split_whitespace().collect();

        match tokens[0].strip_suffix('%') {
            Some(bin) => match bin.parse::<usize>() {
                Ok(b) => target_cnt_by_bin[b] += 1,
                Err(_) => bail!("failed to parse bin from: {}", tokens[0]),
            },
            None => bail!("failed to parse bin from: {}", tokens[0]),
        }
        queries.insert(tokens[1].to_string());
    }

    let queries = queries.iter().collect::<Vec<_>>();
    println!("queries found: {}", queries.len());

    println!("target bin distribution:");
    (BIN_START..=MAX_BIN)
        .collect::<Vec<usize>>()
        .chunks(5)
        .for_each(|bins| {
            bins.iter()
                .for_each(|b| print!("{b}%: {} | ", target_cnt_by_bin[*b]));
            println!();
        });

    let decoy_cnt = (FPR * queries.len() as f32) as usize;

    let mut hits = HashMap::new();
    let name_fn = |p: PathBuf| Some((p.clone(), p.file_stem()?.to_str()?.to_string()));

    for (path, name) in glob::glob(results_dir.join("hmmer*.domtbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.insert(name, hmmer::parse_domtbl(File::open(path)?)?);
    }

    for (path, name) in glob::glob(results_dir.join("mmseqs*.tbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.insert(name, mmseqs::parse_tbl(File::open(path)?)?);
    }

    for (path, name) in glob::glob(results_dir.join("nail*.tbl").to_str().unwrap())?
        .filter_map(Result::ok)
        .filter_map(name_fn)
    {
        hits.insert(name, nail::parse_tbl(File::open(path)?)?);
    }

    let mut plot = Plot::new();

    let layout = Layout::new().template(template());
    plot.set_layout(layout);

    for (name, hit_list) in hits.into_iter() {
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
        let mut hit_list: Vec<Hit> = hit_hash.into_values().collect();

        hit_list.sort_by(|a, b| {
            a.e_value
                .partial_cmp(&b.e_value)
                .expect("NaN E-value encountered")
        });

        let mut num_hits_by_bin = vec![0usize; BIN_START + N_BINS];

        let mut decoys_found = 0;
        for hit in hit_list.iter() {
            match hit.target.contains("decoy") {
                true => {
                    decoys_found += 1;
                    if decoys_found >= decoy_cnt {
                        break;
                    }
                }
                false => {
                    let info = extract_target_info(hit)?;
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

        let (x, y): (Vec<usize>, Vec<f32>) = recall_by_bin
            .into_iter()
            .enumerate()
            .skip(BIN_START)
            .unzip();

        plot.add_trace(Scatter::new(x, y).name(name).mode(Mode::Lines));
    }

    plot.write_image("plot", ImageFormat::SVG, 960, 720, 1.0)
        .expect("Failed to export plot");

    Ok(())
}
