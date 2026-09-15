//! What a planned run of a pipeline will cost.
//!
//! A pipeline is a multiset of searches, and this is the one place that
//! multiset is written down: how many sensitivities, how many cells, and which
//! of them overlap. [`super::fit`] supplies what each one costs; this supplies
//! how many there are. A nail run appears here as its two halves, seeding and
//! then the alignment off those seeds.
//!
//! Only the searches are priced. Splitting the query, building mmseqs'
//! database and reformatting what it found are all real wall clock that no
//! total here includes.
//!
//! The arithmetic follows the same rule the ledger folds by, so the total is
//! wall-clock-shaped rather than a sum of core-seconds: commands that ran at
//! once take the longest of themselves, commands that followed one another
//! add. hmmer is the only batched step here, and its parts run together, so
//! the step costs one part rather than all of them.
//!
//! This is a second description of the pipelines -- `recall.rs`,
//! `cloud_search.rs` and `hit_loss.rs` are the first -- so if one moves
//! without the other the prediction goes wrong with nothing to show for it.
//! The guard is `--against`: hand it a finished pipeline directory and it
//! prints what it predicted beside what that run actually took. A composition
//! that has drifted shows up as an error the fit cannot account for.

use std::path::PathBuf;

use anyhow::{Context, bail};
use clap::Parser;

use util::ledger::{self, Ledger};

use super::fit::Model;
use super::{FILE, Part};
use crate::recall::{MMSEQS_S, NAIL_S};

#[derive(Parser, Debug)]
pub struct Args {
    /// Which pipeline to cost
    #[arg(value_name = "recall|cloud-search|hit-loss")]
    pipeline: String,

    /// How many MGnify sequences one shard holds
    #[arg(long, default_value_t = 1_000_000, value_name = "N")]
    seqs: usize,

    /// Mean residues per sequence. The default is what MGnify runs to
    #[arg(long, default_value_t = 320.0, value_name = "X")]
    seq_len: f64,

    /// How many target shards, for a pipeline with a shard axis
    #[arg(long, default_value_t = 1, value_name = "N")]
    shards: usize,

    /// Threads the planned run will use. It has to be the count the
    /// calibration ran at: the model is fitted against wall clock, which does
    /// not move to another one
    #[arg(short, long, default_value_t = 8, value_name = "N")]
    threads: usize,

    /// The model to cost against
    #[arg(long, value_name = "cost.tbl")]
    cost: Option<PathBuf>,

    /// Check the composition against a pipeline that has actually run: a
    /// directory, or a name under benchmarks/mgy/outputs/
    #[arg(long, value_name = "dir|name")]
    against: Option<String>,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let cost = match &args.cost {
        Some(path) => path.clone(),
        None => crate::outputs().join("calibrate").join(FILE),
    };
    let model = Model::read(&cost)?;

    let plan = Plan {
        pipeline: args.pipeline.clone(),
        t_res: args.seqs as f64 * args.seq_len,
        shards: args.shards,
        jobs: crate::search::jobs(args.threads),
    };

    let steps = plan.steps(&model)?;
    report(&plan, &steps, &args);

    if let Some(name) = &args.against {
        score(&crate::parse::pipeline(name)?, &steps)?;
    }

    Ok(())
}

/// The run being costed.
struct Plan {
    pipeline: String,
    t_res: f64,
    shards: usize,
    /// How many ways hmmer's query is cut.
    //
    // this only widens the memory estimate. every slope in
    // cost.tbl was fitted against wall clock at the thread
    // count the calibration ran at, so the wall clock cannot
    // be moved to another one: predict for the machine you
    // calibrated on, or calibrate again
    jobs: usize,
}

/// One line of the estimate: a primitive, how many of it, and what they come
/// to together.
struct Cost {
    what: String,
    times: usize,
    each_s: f64,
    /// Wall clock for all of them, which is not `each * times` where they
    /// overlap.
    wall_s: f64,
    max_rss_kb: Option<u64>,
    holdout: Option<f64>,
}

impl Plan {
    fn steps(&self, model: &Model) -> anyhow::Result<Vec<Cost>> {
        match self.pipeline.as_str() {
            "recall" => self.recall(model),
            "cloud-search" => self.cloud_search(model),
            "hit-loss" => self.hit_loss(model),
            other => bail!(
                "no composition for {other:?}; it is one of recall, cloud-search, hit-loss"
            ),
        }
    }

    /// One search, `times` of them one after another.
    fn serial(
        &self,
        model: &Model,
        part: Part,
        knobs: &[(&str, &str)],
        times: usize,
        what: &str,
    ) -> anyhow::Result<Cost> {
        let row = model.of(part, knobs).with_context(|| {
            format!("cost.tbl has no fit for {part} at {knobs:?}; was it measured?")
        })?;

        let each = row.at(self.t_res).max(0.0);

        Ok(Cost {
            what: what.to_string(),
            times,
            each_s: each,
            wall_s: each * times as f64,
            max_rss_kb: row.max_rss_kb,
            holdout: row.holdout,
        })
    }

    /// hmmer: the query cut `jobs` ways, every part searched at once.
    fn hmmer(&self, model: &Model, times: usize) -> anyhow::Result<Cost> {
        // the load is the whole comparison, not one part of
        // it: a batched step's ledger row is the longest of
        // its commands, so `fit` regressed against the whole
        // thing and the split is already inside the slope
        let mut cost = self.serial(
            model,
            Part::Hmmer,
            &[],
            times,
            &format!("hmmer ({} parts at once)", self.jobs),
        )?;

        // the parts run together, so the step holds all of their memory
        cost.max_rss_kb = cost.max_rss_kb.map(|rss| rss * self.jobs as u64);
        Ok(cost)
    }

    fn recall(&self, model: &Model) -> anyhow::Result<Vec<Cost>> {
        let mut steps = Vec::new();

        // recall searches every shard once per sensitivity, and the two tools
        // are swept over different ones, so a single nail row and a single
        // mmseqs row would price a fifth of what runs.
        //
        // a nail run is its seeding plus the alignment off those seeds; the
        // sweep times the halves rather than the whole, since that is what
        // cloud-search and hit-loss are built out of
        for s in NAIL_S.split(',') {
            for part in [Part::Seed, Part::Align] {
                steps.push(self.serial(
                    model,
                    part,
                    &[("s", s)],
                    self.shards,
                    &format!("nail {part} -s {s} (per shard)"),
                )?);
            }
        }

        for s in MMSEQS_S.split(',') {
            steps.push(self.serial(
                model,
                Part::Mmseqs,
                &[("s", s)],
                self.shards,
                &format!("mmseqs -s {s} (per shard)"),
            )?);
        }

        steps.push(self.hmmer(model, self.shards)?);
        Ok(steps)
    }

    fn cloud_search(&self, model: &Model) -> anyhow::Result<Vec<Cost>> {
        // the default grid: nine alphas by nine betas, and the unpruned cell
        const CELLS: usize = 9 * 9 + 1;

        // every cell replays the one seed set, so the alignment is priced at
        // the seeding sensitivity cloud-search uses rather than at recall's
        let seed_s = crate::cloud_search::MMSEQS_S;

        Ok(vec![
            self.serial(
                model,
                Part::Seed,
                &[("s", seed_s)],
                1,
                "seed (once, replayed by every cell)",
            )?,
            self.hmmer(model, 1)?,
            self.serial(
                model,
                Part::Align,
                &[("s", seed_s)],
                CELLS,
                &format!("align ({CELLS} cells)"),
            )?,
        ])
    }

    fn hit_loss(&self, model: &Model) -> anyhow::Result<Vec<Cost>> {
        let seed_s = crate::hit_loss::MMSEQS_S;

        Ok(vec![
            self.serial(model, Part::Seed, &[("s", seed_s)], 1, "seed")?,
            self.hmmer(model, 1)?,
            self.serial(model, Part::Align, &[("s", seed_s)], 1, "align")?,
        ])
    }
}

fn report(plan: &Plan, steps: &[Cost], args: &Args) {
    let total: f64 = steps.iter().map(|step| step.wall_s).sum();
    let peak = steps.iter().filter_map(|step| step.max_rss_kb).max();

    println!(
        "{} over {} shard(s) of {} sequences, against all of Pfam, at {} threads\n",
        plan.pipeline, plan.shards, args.seqs, args.threads
    );

    println!("  {:<44} {:>5}  {:>12}  {:>12}", "step", "n", "each", "wall");
    println!("  {:-<44} {:->5}  {:->12}  {:->12}", "", "", "", "");

    for step in steps {
        println!(
            "  {:<44} {:>5}  {:>12}  {:>12}",
            step.what,
            step.times,
            clock(step.each_s),
            clock(step.wall_s),
        );
    }

    println!("  {:-<44} {:->5}  {:->12}  {:->12}", "", "", "", "");
    println!("  {:<44} {:>5}  {:>12}  {:>12}", "total", "", "", clock(total));

    if let Some(peak) = peak {
        println!("\n  peak resident set, one step at a time: {:.1} GB", peak as f64 / (1 << 20) as f64);
    }

    // the worst holdout is what the whole estimate is worth: a model that
    // missed the rung it did not see will miss this by at least as much
    if let Some(worst) = steps
        .iter()
        .filter_map(|step| step.holdout)
        .max_by(|x, y| x.abs().total_cmp(&y.abs()))
    {
        println!(
            "  the worst-extrapolating part was {worst:+.1}% out on the rung it \
             was not fitted on"
        );
    }
}

/// What the composition predicted, beside what a real run took.
fn score(dir: &std::path::Path, steps: &[Cost]) -> anyhow::Result<()> {
    let ran = Ledger::load(dir)?;
    ledger::warn(ran.failed(), "command(s)");

    let actual: f64 = ran.rows().filter_map(|row| row.wall_s).sum();
    let predicted: f64 = steps.iter().map(|step| step.wall_s).sum();

    println!("\n  against {}", dir.display());
    println!("    predicted {:>12}", clock(predicted));
    println!("    actual    {:>12}", clock(actual));

    if actual > 0.0 {
        println!(
            "    out by    {:>11.1}%",
            100.0 * (predicted - actual) / actual
        );
    }

    Ok(())
}

/// Seconds, in whatever unit makes the number readable.
fn clock(seconds: f64) -> String {
    match seconds {
        s if s < 90.0 => format!("{s:.1}s"),
        s if s < 5400.0 => format!("{:.1}m", s / 60.0),
        s if s < 172_800.0 => format!("{:.1}h", s / 3600.0),
        s => format!("{:.1}d", s / 86_400.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_clock_reads_in_the_unit_that_suits_it() {
        assert_eq!(clock(12.0), "12.0s");
        assert_eq!(clock(600.0), "10.0m");
        assert_eq!(clock(7200.0), "2.0h");
        assert_eq!(clock(432_000.0), "5.0d");
    }
}
