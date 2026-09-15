//! Turning a calibrate ledger into a cost model.
//!
//! One row of `cost.tbl` per search per parameterization: an intercept, a
//! slope in target residues, and what the fit is worth.
//!
//! ```text
//! seconds = intercept + slope * t_res
//! ```
//!
//! The intercept is what the command pays before it starts, and separating it
//! is not pedantry. Fitting a proportion -- seconds per residue, averaged --
//! folds the fixed cost into a slope and then multiplies it by the size of the
//! run, which at a million targets turns a one-minute startup into an hour
//! that is not there. nail has a large one, since it builds an mmseqs profile
//! database out of the HMMs on every invocation and the query never changes.
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

use super::run::PART;
use super::{FILE, Part};
use crate::inputs;

/// The one term a cost is fitted against, and its column in `cost.tbl`.
const TERM: &str = "t_res";

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

    let sizes = inputs::ladder::target_residues()?;
    let points = points(&ran, &sizes)?;
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
    // s, A, B -- an alignment at -A 2 -B 4 is a different cost
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
    t_res: f64,
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

/// Every ledger row that is a timed search, turned into a measurement.
fn points(ran: &Ledger, sizes: &BTreeMap<usize, u64>) -> anyhow::Result<Vec<Point>> {
    let mut out = Vec::new();

    for row in ran.rows() {
        let Some(part) = row.params.get(PART).and_then(|key| Part::parse(key)) else {
            continue;
        };
        let Some(wall_s) = row.wall_s else { continue };

        let t_rung: usize = row
            .shard
            .parse()
            .with_context(|| format!("shard {:?} is not a target rung", row.shard))?;

        let t_res = *sizes
            .get(&t_rung)
            .with_context(|| format!("no target rung {t_rung} in sizes.tbl"))?
            as f64;

        out.push(Point {
            key: Key {
                part,
                knobs: knobs(&row.params),
            },
            t_rung,
            t_res,
            wall_s,
            cpu_s: row.cpu_s,
            max_rss_kb: row.max_rss_kb,
        });
    }

    Ok(out)
}

/// Everything that tells one measurement of a primitive from another, which is
/// every param but the two that say where on the ladder it sat.
fn knobs(params: &BTreeMap<String, String>) -> BTreeMap<String, String> {
    params
        .iter()
        .filter(|(key, _)| key.as_str() != PART)
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect()
}

fn group(points: &[Point]) -> BTreeMap<Key, Vec<Point>> {
    let mut out: BTreeMap<Key, Vec<Point>> = BTreeMap::new();
    for point in points {
        out.entry(point.key.clone()).or_default().push(point.clone());
    }
    out
}

/// How many genuinely different sizes the points cover.
fn distinct(points: &[Point]) -> usize {
    let mut sizes: Vec<u64> = points.iter().map(|p| p.t_res.to_bits()).collect();

    sizes.sort_unstable();
    sizes.dedup();
    sizes.len()
}

/// One row of the model.
struct Fit {
    key: Key,
    /// The intercept, then the slope in target residues.
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
        format!("{:.3} + {:.3e}*{TERM}", self.coef[0], self.coef[1])
    }

    fn of(key: &Key, points: &[Point], held: Option<usize>) -> Option<Fit> {
        // two rungs at least, or there is no line: infinitely many pass
        // through one point, and points that all sat at one size fit a slope
        // of anything
        if distinct(points) < 2 {
            return None;
        }

        let all = line(points, |p| p.wall_s)?;

        // fit again without the top rung, then evaluate that fit on the
        // rung it excluded
        let holdout = held.and_then(|held| {
            let (kept, out): (Vec<Point>, Vec<Point>) =
                points.iter().cloned().partition(|p| p.t_rung != held);

            if out.is_empty() || distinct(&kept) < 2 {
                return None;
            }

            let short = line(&kept, |p| p.wall_s)?;
            let want: f64 = out.iter().map(|p| p.wall_s).sum::<f64>() / out.len() as f64;
            let got: f64 = out.iter().map(|p| short.at(p.t_res)).sum::<f64>() / out.len() as f64;

            (want > 0.0).then(|| 100.0 * (got - want) / want)
        });

        Some(Fit {
            key: key.clone(),
            coef: vec![all.intercept, all.slope],
            r2: all.r2,
            holdout,
            max_rss_kb: points.iter().filter_map(|p| p.max_rss_kb).max(),
            n: points.len(),
        })
    }
}

/// A least-squares fit of `y = intercept + slope * t_res`.
struct Line {
    intercept: f64,
    slope: f64,
    r2: f64,
}

impl Line {
    fn at(&self, t_res: f64) -> f64 {
        self.intercept + self.slope * t_res
    }
}

/// Least squares of the points' cost against their target residues.
fn line(points: &[Point], y_of: impl Fn(&Point) -> f64) -> Option<Line> {
    let n = points.len() as f64;
    if points.len() < 2 {
        return None;
    }

    let mean_x = points.iter().map(|p| p.t_res).sum::<f64>() / n;
    let mean_y = points.iter().map(&y_of).sum::<f64>() / n;

    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for point in points {
        let dx = point.t_res - mean_x;
        let dy = y_of(point) - mean_y;

        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }

    // every rung at one size: no spread to take a slope from
    if sxx < 1e-12 {
        return None;
    }

    let slope = sxy / sxx;

    Some(Line {
        intercept: mean_y - slope * mean_x,
        slope,
        // a flat measurement has no spread to account for, and is a perfect
        // fit by convention
        r2: match syy == 0.0 {
            true => 1.0,
            false => (sxy * sxy) / (sxx * syy),
        },
    })
}

fn write(out: &Path, fits: &[Fit], held: Option<usize>) -> anyhow::Result<()> {
    let mut headers = vec!["part".to_string(), "knobs".to_string(), "intercept".to_string()];
    headers.push(TERM.to_string());
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

            let mut cells = vec![
                fit.key.part.to_string(),
                match knobs.is_empty() {
                    true => "-".to_string(),
                    false => knobs.join(","),
                },
                format!("{:.4}", fit.coef[0]),
                format!("{:.6e}", fit.coef[1]),
            ];

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
        "# seconds = intercept + slope * t_res\n\
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
    /// The intercept, then the slope in target residues.
    pub coef: Vec<f64>,
    pub max_rss_kb: Option<u64>,
    pub holdout: Option<f64>,
}

impl Row {
    /// What this search costs against a target of this many residues.
    pub fn at(&self, t_res: f64) -> f64 {
        self.coef[0] + self.coef[1] * t_res
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

            let coef = vec![
                cell("intercept").parse().unwrap_or(0.0),
                cell(TERM).parse().unwrap_or(0.0),
            ];

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

    fn point(t_res: f64, wall_s: f64, t_rung: usize) -> Point {
        Point {
            key: Key {
                part: Part::Align,
                knobs: BTreeMap::new(),
            },
            t_rung,
            t_res,
            wall_s,
            cpu_s: None,
            max_rss_kb: None,
        }
    }

    /// The point of separating the intercept: nail builds an mmseqs profile
    /// database out of the HMMs before it looks at a target, so a fit that
    /// smeared that into the slope would multiply it by the size of the run.
    #[test]
    fn a_line_recovers_an_intercept_and_a_slope() {
        let points: Vec<Point> = (1..=4)
            .map(|i| point(i as f64, 10.0 + 2.0 * i as f64, i))
            .collect();

        let fit = line(&points, |p| p.wall_s).unwrap();

        assert!((fit.intercept - 10.0).abs() < 1e-9, "{}", fit.intercept);
        assert!((fit.slope - 2.0).abs() < 1e-9, "{}", fit.slope);
        assert!((fit.r2 - 1.0).abs() < 1e-9);
    }

    /// One size, however many times it was measured, says nothing about how
    /// the cost grows -- so it is refused rather than fitted to a slope of
    /// whatever the noise happened to be.
    #[test]
    fn one_size_is_not_enough_to_fit() {
        let points = vec![point(5.0, 1.0, 1), point(5.0, 1.2, 1)];

        assert!(line(&points, |p| p.wall_s).is_none());
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

    /// A cost that bends away from the line is caught by the holdout, which is
    /// the failure the rungs were made to double for.
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
