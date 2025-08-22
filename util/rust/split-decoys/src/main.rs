use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

use bioio::fasta::Fasta;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        println!("usage: split-decoys <profmark.test.fa>");
        return Ok(());
    };

    let fa_path = std::path::Path::new(&args[1]);
    let mut positives_out = BufWriter::new(File::create(fa_path.with_file_name("positives.fa"))?);
    let mut decoys_out = BufWriter::new(File::create(fa_path.with_file_name("decoys.fa"))?);

    let (positives, decoys): (Vec<_>, Vec<_>) = Fasta::from_path(fa_path)?
        .records
        .into_values()
        .partition(|r| !r.name.starts_with("decoy"));

    positives
        .iter()
        .try_for_each(|r| writeln!(positives_out, "{r}"))?;

    decoys
        .iter()
        .try_for_each(|r| writeln!(decoys_out, "{r}"))?;

    Ok(())
}
