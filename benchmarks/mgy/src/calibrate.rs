//! What a search costs here, and what a planned run will.
//!
//! The other pipelines ask what was found. This one asks what it took, and it
//! is the only pipeline whose artifact is a model rather than a table of hits.
//! A run of Pfam against a million MGnify sequences takes long enough that
//! finding out by running it is not an answer, so this measures the searches
//! at sizes that fit in an afternoon and composes them.
//!
//! ```text
//! calibrate run      time every search over the rungs of the ladder
//! calibrate fit      turn those timings into cost.tbl
//! calibrate predict  compose a pipeline out of cost.tbl
//! ```
//!
//! The query is every Pfam family and stays that way, because that is what the
//! benchmarks search with. So the only thing that moves is the target, and a
//! part's cost is a line in target residues:
//!
//! ```text
//! seconds = intercept + slope * t_res
//! ```
//!
//! An earlier version of this swept a grid on both axes and fitted a bilinear
//! surface. Holding the query fixed collapses three of those four terms into
//! the two here, and the terms it drops were the ones carrying the error: the
//! product coefficient did most of the work and moved by a factor of two when
//! a single rung left the grid.
//!
//! What makes this a calibration rather than a stopwatch is that [`fit`] holds
//! the top rung back, fits on the rest, and scores its own prediction of the
//! rung it did not see, so a model that extrapolates badly reports it in
//! `cost.tbl` instead of being believed.
//!
//! `mgy cutoffs` is the other thing this crate calls a calibration, and the
//! word does not mean the same: cutoffs calibrates scores, what a hit has to
//! beat to count, and this calibrates cost.

pub mod fit;
pub mod predict;
pub mod run;

use clap::{Subcommand, ValueEnum};

/// The file a fitted cost model is written to.
pub const FILE: &str = "cost.tbl";

#[derive(Subcommand)]
pub enum Cmd {
    /// Time every search over the rungs of the ladder.
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

/// The searches a pipeline is built out of.
///
/// Only searches. Building a database, cutting the query up and reformatting
/// what came back are all real wall clock, but they are not what the
/// benchmarks compare, and pricing them here made the model answer a question
/// nobody asked.
///
/// A whole nail search is [`Part::Seed`] then [`Part::Align`], so it is not a
/// part of its own. Measured over a ladder spanning 131x, the two add up to
/// the end-to-end run within 1.1% at every sensitivity, which is what retired
/// the third timing.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Part {
    /// `nail --only-seed`: the pass that produces pairs rather than scores.
    Seed,
    /// `nail --seeds`: cloud search and alignment over a seed set that already
    /// exists.
    Align,
    /// `mmseqs search` against a prebuilt database.
    Mmseqs,
    /// `hmmsearch` over every part of the split query at once.
    Hmmer,
}

impl Part {
    pub const ALL: [Part; 4] = [Part::Seed, Part::Align, Part::Mmseqs, Part::Hmmer];

    /// How this part is spelled in `cost.tbl` and in a calibrate ledger.
    pub fn key(self) -> &'static str {
        match self {
            Part::Seed => "seed",
            Part::Align => "align",
            Part::Mmseqs => "mmseqs",
            Part::Hmmer => "hmmer",
        }
    }

    pub fn parse(key: &str) -> Option<Part> {
        Part::ALL.into_iter().find(|part| part.key() == key)
    }
}

impl std::fmt::Display for Part {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.key())
    }
}
