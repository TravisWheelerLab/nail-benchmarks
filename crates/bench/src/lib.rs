//! What the benchmarks in this repo agree on.
//!
//! Each benchmark asks its own question and keeps its own analyses. What they
//! share is the shape of the record a run leaves behind: a [`pail::Table`]
//! sink writes `manifest.tbl`, `parse` reads it back through [`manifest`], and
//! whatever it works out gets written through [`tbl`].
//!
//! Nothing here knows what a hit is or what counts as a true one. That is the
//! part every benchmark answers differently, and it stays with the benchmark.

pub mod manifest;
pub mod tbl;
