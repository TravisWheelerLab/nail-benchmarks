//! What a search costs here, and what a planned run will.
//!
//! The other pipelines ask what was found. This one asks what it took, and it
//! is the only pipeline whose artifact is a model rather than a table of hits.
//! A run of Pfam against a million MGnify sequences takes long enough that
//! finding out by running it is not an answer, so this measures the parts at
//! sizes that fit in an afternoon and composes them.
//!
//! ```text
//! calibrate run      time every primitive over the ladder
//! calibrate fit      turn those timings into cost.tbl
//! calibrate predict  compose a pipeline out of cost.tbl
//! ```
//!
//! What makes this a calibration rather than a stopwatch is that [`fit`] holds
//! the top rung back, fits on the rest, and scores its own prediction of the
//! rung it did not see, so a model that extrapolates badly reports it in
//! `cost.tbl` instead of being believed.
//!
//! A pipeline is built out of a handful of distinct commands, and every
//! pipeline here is some multiset of them. [`run`] times each one using the
//! same builders in [`crate::search`] that the real pipelines use, so what is
//! measured is what will run rather than a model of it. The replay family is
//! the one that does not follow the others: cloud-search seeds once and
//! replays every `(A, B)` cell off that one seed set, so a cell's cost follows
//! the seed count rather than the size of the target. That makes the seed
//! count a bridge, which [`fit`] learns as a function of the search size and
//! feeds to the replay surface. Without it there is no predicting a cell at a
//! size nobody has seeded.
//!
//! `mgy cutoffs` is the other thing this crate calls a calibration, and the
//! word does not mean the same: cutoffs calibrates scores, what a hit has to
//! beat to count, and this calibrates cost.

pub mod fit;
pub mod predict;
pub mod run;

use clap::Subcommand;

/// The file a fitted cost model is written to.
pub const FILE: &str = "cost.tbl";

#[derive(Subcommand)]
pub enum Cmd {
    /// Time every primitive over the rungs of the ladder.
    Run(run::Args),
    /// Fit the timings into a cost model, and score it against a held-out
    /// rung.
    Fit(fit::Args),
    /// What a planned run of a pipeline will cost.
    Predict(predict::Args),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Run(args) => run::main(args),
        Cmd::Fit(args) => fit::main(args),
        Cmd::Predict(args) => predict::main(args),
    }
}

/// The parts a pipeline is built out of: one command, or one batch of them,
/// that some pipeline runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Part {
    /// Cutting the query set into parts for hmmer.
    Split,
    /// Building mmseqs' database over a target.
    Createdb,
    /// `nail --only-seed`: the pass that produces pairs rather than scores.
    Seed,
    /// `nail --seeds`: cloud search and alignment over a seed set that already
    /// exists. cloud-search's cells and hit-loss's one run are both this.
    Replay,
    /// `nail search` end to end, seeding inside itself. recall's column.
    Nail,
    /// `mmseqs search` against a prebuilt database.
    Mmseqs,
    /// `mmseqs convertalis`, which reformats what the search found.
    Convert,
    /// Not a command: how many pairs seeding found.
    //
    // a replay's cost follows the seed set it replays rather
    // than the target that produced it, so predicting a cell
    // at a size nobody has seeded means predicting the
    // seeding first. fitted like any other part, against the
    // same load; its slope is pairs rather than seconds
    Seeds,
    /// One `hmmsearch` over one part of the split query, plus the `cat` that
    /// gathers the parts back up.
    Hmmer,
}

impl Part {
    pub const ALL: [Part; 9] = [
        Part::Split,
        Part::Createdb,
        Part::Seed,
        Part::Replay,
        Part::Nail,
        Part::Mmseqs,
        Part::Convert,
        Part::Seeds,
        Part::Hmmer,
    ];

    /// How this part is spelled in `cost.tbl` and in a calibrate ledger.
    pub fn key(self) -> &'static str {
        match self {
            Part::Split => "split",
            Part::Createdb => "createdb",
            Part::Seed => "seed",
            Part::Replay => "replay",
            Part::Nail => "nail",
            Part::Mmseqs => "mmseqs",
            Part::Convert => "convert",
            Part::Seeds => "seeds",
            Part::Hmmer => "hmmer",
        }
    }

    pub fn parse(key: &str) -> Option<Part> {
        Part::ALL.into_iter().find(|part| part.key() == key)
    }

    /// The terms this part's cost is fitted against.
    pub fn loads(self) -> &'static [Load] {
        match self {
            // cutting a file up is proportional to the file
            Part::Split => &[Load::Query],
            // so is reading one into a database
            Part::Createdb => &[Load::Target],

            // a search needs all three, which is why the ladder sweeps a
            // grid: the product alone is not a stand-in for its cost, at
            // least not for nail and mmseqs. both pay for the query before
            // they look at a target, since nail builds a profile database
            // out of the HMMs every run, and both pay for the target
            // whatever the query is -- an index to build, a database to
            // scan. only what is left over is the comparison
            //
            // measured over the target axis, at two query rungs 4.69x
            // apart, the slope in target residues rose:
            //     nail    1.93x
            //     mmseqs  2.59x
            //     hmmer   3.75x
            // a pure product would have moved all three by 4.69x, and a
            // query-independent cost by 1.0. fitting only query*target
            // forces the other two terms into it, and the error is then
            // multiplied by the size of the run
            Part::Seed | Part::Nail | Part::Mmseqs | Part::Hmmer | Part::Convert => {
                &[Load::Query, Load::Target, Load::Product]
            }

            // how many pairs come back is mostly the comparison, but a
            // bigger target offers more to match against at any query size
            Part::Seeds => &[Load::Query, Load::Target, Load::Product],

            // the target reaches a replay through the seed count rather
            // than directly, so a term of its own would restate one it
            // already has
            Part::Replay => &[Load::Query, Load::Seeds],
        }
    }
}

impl std::fmt::Display for Part {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.key())
    }
}

/// What a part's cost is proportional to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Load {
    /// Query residues.
    //
    // residues rather than sequence counts throughout:
    // MGnify's lengths vary enough that a count is the wrong
    // x-axis, and build already writes the residues down in
    // sizes.tbl
    Query,
    /// Target residues.
    Target,
    /// Query residues times target residues: the cells of the comparison.
    Product,
    /// How many pairs seeding found.
    Seeds,
}

impl Load {
    /// Every term any part can have, in the order `cost.tbl` columns them.
    pub const ALL: [Load; 4] = [Load::Query, Load::Target, Load::Product, Load::Seeds];

    pub fn key(self) -> &'static str {
        match self {
            Load::Query => "q_res",
            Load::Target => "t_res",
            Load::Product => "q_res*t_res",
            Load::Seeds => "seeds",
        }
    }

    /// The load itself, given the two axes and what a search off them found.
    pub fn of(self, q_res: f64, t_res: f64, seeds: f64) -> f64 {
        match self {
            Load::Query => q_res,
            Load::Target => t_res,
            Load::Product => q_res * t_res,
            Load::Seeds => seeds,
        }
    }
}
