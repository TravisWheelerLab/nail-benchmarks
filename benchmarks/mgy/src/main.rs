//! Pfam against MGnify, three ways.
//!
//! One `build` draws the query set and the target shards; every pipeline
//! searches those same files. What a pipeline owns is everything downstream of
//! them -- its own seeds, its own hmmer truth, its own results -- so a pipeline
//! directory can be read on its own without asking what else has been run.

mod build;
mod cloud_search;
mod cutoffs;
mod hit_loss;
mod manifest;
mod parse;
mod recall;
mod scores;
mod search;

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
    /// Draw the query set and the target shards every pipeline searches.
    Build(build::Args),
    /// Learn per-family false-positive score cutoffs from reversed decoys.
    #[command(subcommand)]
    Cutoffs(cutoffs::Cmd),
    /// Search nail and mmseqs against every shard, sweeping their prefilter
    /// sensitivity.
    Recall(recall::Args),
    /// Seed once, then search every (A, B) cell off those seeds.
    CloudSearch(cloud_search::Args),
    /// Seed once, run hmmer for truth, then run nail once at its defaults.
    HitLoss(hit_loss::Args),
    /// Turn any finished pipeline into its scores table.
    Parse(parse::Args),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Cutoffs(cmd) => cutoffs::main(cmd),
        Command::Recall(args) => recall::main(args),
        Command::CloudSearch(args) => cloud_search::main(args),
        Command::HitLoss(args) => hit_loss::main(args),
        Command::Parse(args) => parse::main(args),
    }
}

/// This crate's directory, fixed at compile time. The shared inputs and every
/// pipeline's output hang off it.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn queries() -> PathBuf {
    dir().join("queries")
}

pub fn targets() -> PathBuf {
    dir().join("targets")
}
