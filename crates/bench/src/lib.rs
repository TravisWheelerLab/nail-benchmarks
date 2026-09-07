//! What the benchmarks in this repo agree on.
//!
//! Each benchmark asks its own question and keeps its own analyses. Two things
//! are the same for all of them. The first is the shape of the record a run
//! leaves behind: a [`pail::Table`] sink writes `manifest.tbl`, `parse` reads it
//! back through [`manifest`], and whatever it works out gets written through
//! [`tbl`]. The second is where the programs being benchmarked and the sequence
//! data they run on were put, which is [`tools`].
//!
//! Nothing here knows what a hit is or what counts as a true one. That is the
//! part every benchmark answers differently, and it stays with the benchmark.

pub mod manifest;
pub mod tbl;
pub mod tools;
