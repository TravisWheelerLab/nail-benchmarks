//! Turning a calibrate ledger into a cost model.
//!
//! One row of `cost.tbl` per primitive per parameterization: an intercept, a
//! coefficient for each term of the part's form, and what the fit is worth.
//!
//! ```text
//! seconds = intercept + a*q_res + b*t_res + c*q_res*t_res
//! ```
//!
//! The intercept is what the command pays before it starts, and separating it
//! is not pedantry. Fitting a proportion -- seconds per residue, averaged --
//! folds the fixed cost into a slope and then multiplies it by the size of the
//! run, which at a million targets turns a one-minute startup into an hour
//! that is not there. Free parameters and a ladder wide enough to fit them on
//! is what the nested rungs are for.
//!
//! Every fit here is used by extrapolating past the data that made it, which
//! is the one thing a regression cannot promise. So the top rung is held back,
//! the fit is made on the rest, and the model is scored on the rung it did not
//! see. `cost.tbl` carries that error per primitive, and it is the number to
//! read before believing a prediction.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::Parser;

use util::ledger::{self, Ledger};
use util::tbl;

use super::run::{PART, QUERY, SEEDS};
use super::{FILE, Load, Part};
use crate::inputs;

#[derive(Parser, Debug)]
pub struct Args {
    /// The calibrate directory, or the name of one under
    /// benchmarks/mgy/outputs/
    #[arg(default_value = "calibrate", value_name = "dir|name")]
    pipeline: String,

    /// Where the model goes. Defaults to cost.tbl beside the ledger
    #[arg(short, long, value_name = "cost.tbl")]
    out: Option<PathBuf>,

    /// Fit on every rung, scoring nothing. The default holds the largest
    /// target rung back and reports the error on it
    #[arg(long)]
    no_holdout: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let dir = crate::parse::pipeline(&args.pipeline)?;
    let out = args.out.unwrap_or_else(|| dir.join(FILE));

    let ran = Ledger::load(&dir)?;
    ledger::warn(ran.failed(), "command(s)");

    let sizes = Sizes::read()?;
    let points = points(&ran, &sizes, &dir.join("results"))?;
    ensure!(!points.is_empty(), "no timed primitives in {}", dir.display());

    let held = match args.no_holdout {
        true => None,
        false => points.iter().map(|point| point.t_rung).max(),
    };

    let mut fits: Vec<Fit> = Vec::new();
    for (key, group) in group(&points) {
        match Fit::of(&key, &group, held) {
            Some(fit) => fits.push(fit),
            None => eprintln!(
                "warning: {} has {} measurement(s) at {} distinct load(s); \
                 it needs two to separate a slope from an intercept",
                key.label(),
                group.len(),
                distinct(&group),
            ),
        }
    }

    ensure!(!fits.is_empty(), "nothing had enough rungs to fit");
    write(&out, &fits, held)?;

    println!("wrote {}", out.display());
    for fit in &fits {
        println!(
            "  {:<28} {:<46} r2={:.4}{}",
            fit.key.label(),
            fit.form(),
            fit.r2,
            match fit.holdout {
                Some(err) => format!("  holdout {err:+.1}%"),
                None => String::new(),
            }
        );
    }

    Ok(())
}

/// One primitive at one parameterization: what `cost.tbl` has a row for.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Key {
    part: Part,
    /// The knobs that tell one measurement of this primitive from another.
    //
    // s, A, B -- a replay at -A 2 -B 4 is a different cost
    // from one at -A 40 -B 64, and one model cannot carry
    // both
    knobs: BTreeMap<String, String>,
}

impl Key {
    fn label(&self) -> String {
        let knobs: String = self
            .knobs
            .iter()
            .map(|(k, v)| format!(" {k}={v}"))
            .collect();

        format!("{}{knobs}", self.part)
    }
}

/// One measurement: what it cost, and how much work it was.
#[derive(Clone, Debug)]
struct Point {
    key: Key,
    /// The target rung, which is what the holdout is taken over.
    t_rung: usize,
    /// One per term of the part's form, in [`Part::loads`] order.
    loads: Vec<f64>,
    /// Seconds, for every part but [`Part::Seeds`], where it is a count of
    /// pairs. The fit is the same line either way.
    wall_s: f64,
    /// Core-seconds.
    //
    // not fitted -- predict answers in wall clock, at the
    // thread count the calibration ran at -- but carried so a
    // reader of ledger.tbl can see the two beside each other
    #[allow(dead_code)]
    cpu_s: Option<f64>,
    max_rss_kb: Option<u64>,
}

/// The residues behind every rung of the ladder, off what `build` wrote.
struct Sizes {
    queries: BTreeMap<usize, u64>,
    targets: BTreeMap<usize, u64>,
}

impl Sizes {
    fn read() -> anyhow::Result<Sizes> {
        Ok(Sizes {
            queries: rungs(&inputs::ladder::sizes(&inputs::ladder::queries()))?,
            targets: rungs(&inputs::ladder::sizes(&inputs::ladder::targets()))?,
        })
    }
}

/// `rung residues bytes`, as `build ladder` writes it.
fn rungs(path: &Path) -> anyhow::Result<BTreeMap<usize, u64>> {
    let table = tbl::read(path)
        .with_context(|| format!("failed to read {}; has `mgy build ladder` run?", path.display()))?;

    let mut out = BTreeMap::new();
    for cells in &table.cells {
        let get = |key: &str| -> anyhow::Result<u64> {
            cells
                .get(key)
                .with_context(|| format!("{} has no {key} column", path.display()))?
                .parse()
                .with_context(|| format!("{} has a {key} that is not a number", path.display()))
        };

        out.insert(get("rung")? as usize, get("residues")?);
    }

    ensure!(!out.is_empty(), "no rungs in {}", path.display());
    Ok(out)
}

/// Every ledger row that is a primitive, turned into a measurement.
fn points(ran: &Ledger, sizes: &Sizes, results: &Path) -> anyhow::Result<Vec<Point>> {
    let mut out = Vec::new();

    for row in ran.rows() {
        // a run says which primitive it is; a stage is named by being one,
        // and one that runs several times per shard suffixes that name so its
        // rows do not collide
        let part = match row.stage.is_empty() {
            true => row.params.get(PART).and_then(|key| Part::parse(key)),
            false => Part::parse(row.stage.split('.').next().unwrap_or(&row.stage)),
        };

        let Some(part) = part else { continue };
        let Some(wall_s) = row.wall_s else { continue };

        // the split is per query rung and searches nothing, so it has no shard
        let t_rung: usize = match row.shard.is_empty() {
            true => 0,
            false => row
                .shard
                .parse()
                .with_context(|| format!("shard {:?} is not a target rung", row.shard))?,
        };

        // a part whose cost has no query term never names one -- building
        // mmseqs' target database is independent of what will be searched
        // against it -- so a missing rung is only a fault where the form
        // needs it
        let wants_query = part
            .loads()
            .iter()
            .any(|load| matches!(load, Load::Query | Load::Product));

        let q_res = match row.params.get(QUERY) {
            Some(q) => {
                let rung: usize = q
                    .parse()
                    .with_context(|| format!("{QUERY}={q:?} is not a query rung"))?;

                *sizes
                    .queries
                    .get(&rung)
                    .with_context(|| format!("no query rung {rung} in sizes.tbl"))?
                    as f64
            }
            None if !wants_query => 0.0,
            None => continue,
        };

        let t_res = match t_rung {
            0 => 0.0,
            rung => *sizes
                .targets
                .get(&rung)
                .with_context(|| format!("no target rung {rung} in sizes.tbl"))? as f64,
        };

        let knobs = knobs(&row.params);

        // a replay's cost follows the seed set it replays,
        // which is a property of the seeding that produced it
        // rather than of the size of the target
        let seeds = match part.loads().contains(&Load::Seeds) {
            false => 0.0,
            true => {
                let label = row.params.get(SEEDS).with_context(|| {
                    format!("replay {:?} does not say which seeding it replayed", row.name)
                })?;

                seeds_at(results, label, &row.shard)? as f64
            }
        };

        let loads: Vec<f64> = part
            .loads()
            .iter()
            .map(|load| load.of(q_res, t_res, seeds))
            .collect();

        out.push(Point {
            key: Key {
                part,
                knobs: knobs.clone(),
            },
            t_rung,
            loads,
            wall_s,
            cpu_s: row.cpu_s,
            max_rss_kb: row.max_rss_kb,
        });

        // every seeding is also a measurement of how many pairs it found, on
        // the same axis. that is what lets a replay be placed at a size
        // nobody has seeded
        if part == Part::Seed {
            // `run` labels the stage `seed.<label>` and the list it wrote
            // `seeds.<label>.<shard>`, so the stage says where to find it
            let label = row.stage.split_once('.').map(|(_, rest)| rest);

            if let Some(label) = label {
                out.push(Point {
                    key: Key {
                        part: Part::Seeds,
                        knobs,
                    },
                    t_rung,
                    loads: Part::Seeds
                        .loads()
                        .iter()
                        .map(|load| load.of(q_res, t_res, 0.0))
                        .collect(),
                    wall_s: seeds_at(results, label, &row.shard)? as f64,
                    cpu_s: None,
                    max_rss_kb: None,
                });
            }
        }
    }

    Ok(out)
}

/// Everything that tells one measurement of a primitive from another, which is
/// every param but the two that say where on the ladder it sat.
fn knobs(params: &BTreeMap<String, String>) -> BTreeMap<String, String> {
    params
        .iter()
        .filter(|(key, _)| !matches!(key.as_str(), PART | QUERY | SEEDS))
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect()
}

/// The seed list one seeding left behind, by the label `run` gave it.
fn seeds_at(results: &Path, label: &str, shard: &str) -> anyhow::Result<u64> {
    lines(&results.join(format!("seeds.{label}.{shard}")))
}

fn lines(path: &Path) -> anyhow::Result<u64> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    Ok(text.lines().filter(|line| !line.trim().is_empty()).count() as u64)
}

fn group(points: &[Point]) -> BTreeMap<Key, Vec<Point>> {
    let mut out: BTreeMap<Key, Vec<Point>> = BTreeMap::new();
    for point in points {
        out.entry(point.key.clone()).or_default().push(point.clone());
    }
    out
}

/// How many genuinely different workloads the points cover.
fn distinct(points: &[Point]) -> usize {
    let mut loads: Vec<Vec<u64>> = points
        .iter()
        .map(|p| p.loads.iter().map(|x| x.to_bits()).collect())
        .collect();

    loads.sort_unstable();
    loads.dedup();
    loads.len()
}

/// One row of the model.
struct Fit {
    key: Key,
    /// The intercept, then one coefficient per term of the part's form.
    coef: Vec<f64>,
    r2: f64,
    /// Per cent, signed: how far the fit was out on the rung it did not see.
    holdout: Option<f64>,
    /// The largest resident set anything in this group reached.
    max_rss_kb: Option<u64>,
    n: usize,
}

impl Fit {
    /// The fitted cost as a printable expression.
    fn form(&self) -> String {
        let terms: String = self
            .key
            .part
            .loads()
            .iter()
            .zip(self.coef.iter().skip(1))
            .map(|(load, c)| format!(" + {c:.3e}*{}", load.key()))
            .collect();

        format!("{:.3}{terms}", self.coef[0])
    }

    fn of(key: &Key, points: &[Point], held: Option<usize>) -> Option<Fit> {
        // more workloads than terms, or there is no fit: infinitely many
        // lines pass through one point, and points that all sat at one size
        // fit a slope of anything
        let terms = key.part.loads().len();
        if distinct(points) <= terms {
            return None;
        }

        let all = line(points, |p| p.wall_s)?;

        // fit again without the top rung, then evaluate that fit on the
        // rung it excluded
        let holdout = held.and_then(|held| {
            let (kept, out): (Vec<Point>, Vec<Point>) =
                points.iter().cloned().partition(|p| p.t_rung != held);

            if out.is_empty() || distinct(&kept) <= terms {
                return None;
            }

            let short = line(&kept, |p| p.wall_s)?;
            let want: f64 = out.iter().map(|p| p.wall_s).sum::<f64>() / out.len() as f64;
            let got: f64 =
                out.iter().map(|p| short.at(&p.loads)).sum::<f64>() / out.len() as f64;

            (want > 0.0).then(|| 100.0 * (got - want) / want)
        });

        Some(Fit {
            key: key.clone(),
            coef: all.coef,
            r2: all.r2,
            holdout,
            max_rss_kb: points.iter().filter_map(|p| p.max_rss_kb).max(),
            n: points.len(),
        })
    }
}

/// A least-squares fit of `y = c[0] + c[1]*x1 + ...` over the terms a point
/// carries.
struct Line {
    /// The intercept first, then one coefficient per term.
    coef: Vec<f64>,
    r2: f64,
}

impl Line {
    fn at(&self, loads: &[f64]) -> f64 {
        self.coef
            .iter()
            .skip(1)
            .zip(loads)
            .fold(self.coef[0], |sum, (c, x)| sum + c * x)
    }
}

/// Least squares over an intercept plus however many terms the points carry.
fn line(points: &[Point], y_of: impl Fn(&Point) -> f64) -> Option<Line> {
    let terms = points.first()?.loads.len();
    let width = terms + 1;

    if points.len() < width {
        return None;
    }

    // what: accumulate the normal equations, X'X | X'y
    // why: at most four by four here, so elimination is
    //      enough and nothing needs a decomposition
    //
    // row i of the design matrix is [1, x1, x2 ...]
    let row = |p: &Point| -> Vec<f64> {
        let mut out = Vec::with_capacity(width);
        out.push(1.0);
        out.extend_from_slice(&p.loads);
        out
    };

    let mut a = vec![vec![0.0; width + 1]; width];
    for point in points {
        let x = row(point);
        let y = y_of(point);

        for i in 0..width {
            for j in 0..width {
                a[i][j] += x[i] * x[j];
            }
            a[i][width] += x[i] * y;
        }
    }

    let coef = solve(&mut a)?;

    let n = points.len() as f64;
    let mean_y = points.iter().map(&y_of).sum::<f64>() / n;

    let fit = Line { coef, r2: 0.0 };

    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    for point in points {
        let y = y_of(point);
        ss_res += (y - fit.at(&point.loads)).powi(2);
        ss_tot += (y - mean_y).powi(2);
    }

    Some(Line {
        // a flat measurement has no spread to account for, and is a perfect
        // fit by convention
        r2: match ss_tot == 0.0 {
            true => 1.0,
            false => 1.0 - ss_res / ss_tot,
        },
        ..fit
    })
}

/// Gaussian elimination with partial pivoting, on an augmented matrix.
fn solve(a: &mut [Vec<f64>]) -> Option<Vec<f64>> {
    let n = a.len();

    for col in 0..n {
        let pivot = (col..n).max_by(|x, y| a[*x][col].abs().total_cmp(&a[*y][col].abs()))?;

        // a column with nothing in it means the terms cannot be told apart by
        // these points -- one rung, or two terms that moved together
        if a[pivot][col].abs() < 1e-12 {
            return None;
        }

        a.swap(col, pivot);

        for r in col + 1..n {
            let factor = a[r][col] / a[col][col];
            let (above, below) = a.split_at_mut(r);

            for (cell, pivot) in below[0][col..].iter_mut().zip(&above[col][col..]) {
                *cell -= factor * pivot;
            }
        }
    }

    let mut out = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = a[i][n];
        for j in i + 1..n {
            sum -= a[i][j] * out[j];
        }
        out[i] = sum / a[i][i];
    }

    Some(out)
}

fn write(out: &Path, fits: &[Fit], held: Option<usize>) -> anyhow::Result<()> {
    let mut headers = vec!["part".to_string(), "knobs".to_string(), "intercept".to_string()];
    headers.extend(Load::ALL.iter().map(|load| load.key().to_string()));
    headers.extend(
        ["max_rss", "n", "r2", "holdout"]
            .iter()
            .map(|h| h.to_string()),
    );

    let rows: Vec<Vec<String>> = fits
        .iter()
        .map(|fit| {
            let knobs: Vec<String> = fit
                .key
                .knobs
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect();

            // a term this part does not have reads as a dash rather than a
            // zero: it has no such coefficient, which is not the same as one
            // that came out at nothing
            let loads = fit.key.part.loads();

            let mut cells = vec![
                fit.key.part.to_string(),
                match knobs.is_empty() {
                    true => "-".to_string(),
                    false => knobs.join(","),
                },
                format!("{:.4}", fit.coef[0]),
            ];

            cells.extend(Load::ALL.iter().map(|want| {
                match loads.iter().position(|load| load == want) {
                    Some(at) => format!("{:.6e}", fit.coef[at + 1]),
                    None => "-".to_string(),
                }
            }));

            cells.push(match fit.max_rss_kb {
                Some(rss) => rss.to_string(),
                None => "-".to_string(),
            });
            cells.push(fit.n.to_string());
            cells.push(format!("{:.4}", fit.r2));
            cells.push(match fit.holdout {
                Some(err) => format!("{err:+.2}"),
                None => "-".to_string(),
            });

            cells
        })
        .collect();

    let meta = format!(
        "# seconds = intercept + slope * load\n\
         # holdout {}\n\
         #\n",
        match held {
            Some(rung) => format!(
                "per cent, fitting without target rung {rung} and predicting it"
            ),
            None => "not taken: every rung was fitted".to_string(),
        }
    );

    tbl::write(
        out,
        tbl::Table {
            meta: &meta,
            headers: &headers,
            rows: &rows,
            ragged_last: false,
        },
    )
}

/// The model, read back for [`super::predict`].
pub struct Model {
    pub rows: Vec<Row>,
}

pub struct Row {
    pub part: Part,
    pub knobs: BTreeMap<String, String>,
    /// The intercept, then one coefficient per term of the part's form.
    pub coef: Vec<f64>,
    pub max_rss_kb: Option<u64>,
    pub holdout: Option<f64>,
}

impl Row {
    /// What this part costs at the given workload.
    pub fn at(&self, loads: &[f64]) -> f64 {
        self.coef
            .iter()
            .skip(1)
            .zip(loads)
            .fold(self.coef.first().copied().unwrap_or(0.0), |sum, (c, x)| {
                sum + c * x
            })
    }
}

impl Model {
    pub fn read(path: &Path) -> anyhow::Result<Model> {
        let table = tbl::read(path).with_context(|| {
            format!(
                "failed to read {}; run `mgy calibrate fit` first",
                path.display()
            )
        })?;

        let mut rows = Vec::new();
        for cells in &table.cells {
            let cell = |key: &str| cells.get(key).map(String::as_str).unwrap_or("-");

            let Some(part) = Part::parse(cell("part")) else {
                bail!("{} names a part {:?} nothing knows", path.display(), cell("part"));
            };

            let knobs = match cell("knobs") {
                "-" => BTreeMap::new(),
                text => text
                    .split(',')
                    .filter_map(|pair| pair.split_once('='))
                    .map(|(k, v)| (k.to_string(), v.to_string()))
                    .collect(),
            };

            let mut coef = vec![cell("intercept").parse().unwrap_or(0.0)];
            coef.extend(
                part.loads()
                    .iter()
                    .map(|load| cell(load.key()).parse::<f64>().unwrap_or(0.0)),
            );

            rows.push(Row {
                part,
                knobs,
                coef,
                max_rss_kb: cell("max_rss").parse().ok(),
                holdout: cell("holdout").parse().ok(),
            });
        }

        ensure!(!rows.is_empty(), "no rows in {}", path.display());
        Ok(Model { rows })
    }

    /// The fit for one primitive, preferring one whose knobs match.
    pub fn of(&self, part: Part, knobs: &[(&str, &str)]) -> Option<&Row> {
        let matches = |row: &&Row| {
            knobs
                .iter()
                .all(|(k, v)| row.knobs.get(*k).map(String::as_str) == Some(*v))
        };

        self.rows
            .iter()
            .filter(|row| row.part == part)
            .find(matches)
            .or_else(|| self.rows.iter().find(|row| row.part == part))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A part with one term, so the fixtures stay readable.
    fn point(load: f64, wall_s: f64, t_rung: usize) -> Point {
        Point {
            key: Key {
                part: Part::Createdb,
                knobs: BTreeMap::new(),
            },
            t_rung,
            loads: vec![load],
            wall_s,
            cpu_s: None,
            max_rss_kb: None,
        }
    }

    /// Two terms: a query-proportional setup and the comparison itself, which
    /// is the shape every search here really has.
    fn search(q_res: f64, t_res: f64, wall_s: f64, t_rung: usize) -> Point {
        Point {
            key: Key {
                part: Part::Nail,
                knobs: BTreeMap::new(),
            },
            t_rung,
            loads: vec![q_res, q_res * t_res],
            wall_s,
            cpu_s: None,
            max_rss_kb: None,
        }
    }

    /// The point of separating the terms: a command with a fixed cost and a
    /// proportional one has both recovered, rather than the fixed half being
    /// smeared into the slope.
    #[test]
    fn a_line_recovers_an_intercept_and_a_slope() {
        let points: Vec<Point> = (1..=4)
            .map(|i| point(i as f64, 10.0 + 2.0 * i as f64, i))
            .collect();

        let fit = line(&points, |p| p.wall_s).unwrap();

        assert!((fit.coef[0] - 10.0).abs() < 1e-9, "{:?}", fit.coef);
        assert!((fit.coef[1] - 2.0).abs() < 1e-9, "{:?}", fit.coef);
        assert!((fit.r2 - 1.0).abs() < 1e-9);
    }

    /// A query-proportional setup is recovered as its own term.
    //
    // the failure this form exists to prevent: one scalar
    // intercept pushes the query term into the product, which
    // is then multiplied by the size of the target set
    #[test]
    fn a_query_proportional_setup_does_not_land_in_the_product() {
        let mut points = Vec::new();
        for q in [1.0, 5.0] {
            for (i, t) in [2.0, 5.0, 10.0, 20.0].iter().enumerate() {
                let wall = 0.5 + 3.0 * q + 1e-6 * q * t;
                points.push(search(q, *t, wall, i));
            }
        }

        let fit = line(&points, |p| p.wall_s).unwrap();

        assert!((fit.coef[0] - 0.5).abs() < 1e-6, "{:?}", fit.coef);
        assert!((fit.coef[1] - 3.0).abs() < 1e-6, "{:?}", fit.coef);
        assert!((fit.coef[2] - 1e-6).abs() < 1e-9, "{:?}", fit.coef);
    }

    /// One size, however many times it was measured, says nothing about how
    /// the cost grows -- so it is refused rather than fitted to a slope of
    /// whatever the noise happened to be.
    #[test]
    fn one_load_is_not_enough_to_fit() {
        let points = vec![point(5.0, 1.0, 1), point(5.0, 1.2, 1)];

        assert!(line(&points, |p| p.wall_s).is_none());
        assert!(Fit::of(&points[0].key.clone(), &points, None).is_none());
    }

    /// A ladder too short to separate two terms is refused rather than solved
    /// arbitrarily.
    #[test]
    fn two_terms_need_more_than_two_points() {
        let points = vec![search(1.0, 2.0, 3.5, 0), search(1.0, 5.0, 3.5, 1)];

        assert!(Fit::of(&points[0].key.clone(), &points, None).is_none());
    }

    /// The whole claim of the model: fit without the top rung, then evaluate
    /// that fit on it. A perfectly linear cost scores zero error.
    #[test]
    fn the_holdout_scores_the_rung_that_was_not_fitted() {
        let points: Vec<Point> = (1..=4)
            .map(|i| point(i as f64, 10.0 + 2.0 * i as f64, i))
            .collect();

        let fit = Fit::of(&points[0].key.clone(), &points, Some(4)).unwrap();
        assert!(fit.holdout.unwrap().abs() < 1e-6, "{:?}", fit.holdout);
    }

    /// A cost that bends away from the line is caught by the holdout.
    #[test]
    fn a_holdout_catches_a_model_that_does_not_extrapolate() {
        let mut points: Vec<Point> = (1..=3)
            .map(|i| point(i as f64, 2.0 * i as f64, i))
            .collect();

        // the top rung costs twice what a line through the others predicts
        points.push(point(4.0, 16.0, 4));

        let fit = Fit::of(&points[0].key.clone(), &points, Some(4)).unwrap();
        assert!(fit.holdout.unwrap() < -40.0, "{:?}", fit.holdout);
    }
}
