//! Removing what this benchmark generated.

use clap::Parser;

use crate::inputs;

#[derive(Parser, Debug)]
pub struct Args {
    /// Also remove the profmark split, which costs a create-profmark run over
    /// Pfam and which every rebuild otherwise draws from
    #[arg(long)]
    all: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let mut paths = vec![
        ("inputs", inputs::dir()),
        ("outputs", inputs::outputs()),
        ("tmp", inputs::tmp()),
    ];

    if args.all {
        paths.push(("profmark", inputs::profmark()));
    }

    util::clean::run(&inputs::root(), &paths)
}
