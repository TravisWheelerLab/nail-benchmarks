//! Removing what this benchmark generated.

use clap::Parser;

use crate::{cutoffs, inputs};

#[derive(Parser, Debug)]
pub struct Args {
    /// Also remove the calibrations under cutoffs/, which cost a `mgy
    /// cutoffs` run and are what data/mgy-cutoffs.tbl is promoted from
    #[arg(long)]
    all: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let mut paths = vec![
        ("inputs", inputs::root()),
        ("outputs", crate::outputs()),
        ("tmp", crate::tmp()),
    ];

    if args.all {
        paths.push(("cutoffs", cutoffs::root()));
    }

    util::clean::run(&crate::dir(), &paths)
}
