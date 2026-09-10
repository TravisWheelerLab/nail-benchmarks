//! Removing what this benchmark generated.
//!
//! Which is only what it produced. This benchmark's `inputs/` are symlinks to
//! the sequences checked in under `data/long-seqs/`, so there is nothing
//! there for a clean to take back.

use crate::inputs;

pub fn main() -> anyhow::Result<()> {
    let paths = vec![("outputs", inputs::outputs()), ("tmp", inputs::tmp())];

    util::clean::run(&inputs::dir(), &paths)
}
