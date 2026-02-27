use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use bioio::fasta::FastaByteIndex;

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(value_name = "input.fa")]
    fa_path: PathBuf,

    #[arg(short, long, value_name = "path")]
    out_path: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let start = std::time::Instant::now();
    let args = Args::parse();

    let out_path = if let Some(path) = args.out_path {
        path
    } else {
        args.fa_path.with_extension(".rev.fa")
    };

    let mut out = BufWriter::new(File::create(out_path)?);

    let now = std::time::Instant::now();
    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(File::open(&args.fa_path)?)?;

    println!("indexed {:?} in: {:?}", args.fa_path, now.elapsed());

    for i in 1..=index.size {
        let mut rec = index.get_record(i)?;
        rec.reverse();
        writeln!(out, "{rec}")?;
    }

    println!("reversed {:?} in: {:?}", args.fa_path, start.elapsed());

    Ok(())
}
