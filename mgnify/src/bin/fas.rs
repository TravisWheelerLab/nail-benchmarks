use std::{fs::File, path::PathBuf};

use bioio::fasta::FastaByteIndex;

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(value_name = "input.fa")]
    fa_path: PathBuf,

    #[arg(short, conflicts_with = "p")]
    n: Option<usize>,

    #[arg(short, conflicts_with = "n")]
    p: Option<f32>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let now = std::time::Instant::now();
    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(File::open(args.fa_path)?)?;

    let n_sample = match (args.n, args.p) {
        (Some(n), None) => n,
        (None, Some(p)) => (index.size as f32 * p) as usize,
        _ => panic!(),
    };

    for i in 1..=n_sample {
        let seq = index.get(i)?;
        print!("{seq}");
    }

    Ok(())
}
