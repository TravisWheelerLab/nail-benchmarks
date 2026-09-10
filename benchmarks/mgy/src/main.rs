//! Pfam against MGnify, four ways.
//!
//! `build` cuts the two sources into an input set under `inputs/<kind>/`, and
//! every pipeline that reads a set of that shape searches those same files. What a pipeline owns
//! is everything downstream of them -- its own seeds, its own hmmer run, its
//! own results -- so a directory under `outputs/` can be read on its own
//! without asking what else has been run.
//!
//! There are two shapes, and [`inputs`] is where they are written down.

mod analyze;
mod build;
mod cloud_search;
mod cut;
mod cutoffs;
mod hit_loss;
mod import;
mod inputs;
mod parse;
mod plot;
mod recall;
mod scores;
mod search;
mod search_size;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "mgy",
    about = "pfam against mgnify: recall, cloud search, hit loss"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Cut Pfam and MGnify into an input set for the pipelines to search.
    #[command(subcommand)]
    Build(build::Cmd),
    /// Learn per-family false-positive score cutoffs from reversed decoys.
    #[command(subcommand)]
    Cutoffs(cutoffs::Cmd),
    /// Search nail and mmseqs against every shard, sweeping their prefilter
    /// sensitivity.
    Recall(recall::Args),
    /// Seed once, then search every (A, B) cell off those seeds.
    CloudSearch(cloud_search::Args),
    /// Seed once, run hmmer, then run nail once at its defaults.
    HitLoss(hit_loss::Args),
    /// Time every tool over every rung of the query and target ladders.
    SearchSize(search_size::Args),
    /// Bring in result tables produced somewhere other than here, timed by
    /// the `.time` files that came back with them.
    Import(import::Args),
    /// Turn any finished pipeline into its scores table, and that into
    /// numbers.
    #[command(subcommand)]
    Parse(parse::Cmd),
    /// Draw the pruning heatmaps from summary.tbl.
    Plot(plot::Args),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Build(cmd) => build::main(cmd),
        Command::Cutoffs(cmd) => cutoffs::main(cmd),
        Command::Recall(args) => recall::main(args),
        Command::CloudSearch(args) => cloud_search::main(args),
        Command::HitLoss(args) => hit_loss::main(args),
        Command::SearchSize(args) => search_size::main(args),
        Command::Import(args) => import::main(args),
        Command::Parse(cmd) => parse::main(cmd),
        Command::Plot(args) => plot::main(args),
    }
}

/// This crate's directory, fixed at compile time. The shared inputs and every
/// pipeline's output hang off it.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Where the pipelines write. One directory each, under the one name, so what
/// a pipeline read and what it produced are told apart by which of `inputs/`
/// and `outputs/` a path is under, rather than by remembering the pipeline
/// names.
pub fn outputs() -> PathBuf {
    dir().join("outputs")
}

/// Where every pipeline puts its scratch, one directory each.
///
/// Beside `outputs/` rather than inside it, so what a pipeline produced and
/// what it merely needed on the way are not the same tree: `outputs/<name>/`
/// is the record, and everything under here can be deleted without losing
/// any of it.
pub fn tmp() -> PathBuf {
    dir().join("tmp")
}
