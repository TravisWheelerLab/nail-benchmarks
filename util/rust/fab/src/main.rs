use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use bioio::fasta::FastaByteIndex;
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(value_name = "input.fa")]
    fa_path: PathBuf,

    #[arg(short, long, value_name = "n")]
    n_splits: usize,

    #[arg(short, long, value_name = "n", default_value_t = 67779)]
    seed: usize,

    #[arg(short, long, value_name = "path")]
    out_prefix: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let out_prefix = if let Some(path) = args.out_prefix {
        path
    } else {
        args.fa_path.with_extension("")
    };

    if out_prefix.to_string_lossy().ends_with('/') {
        std::fs::create_dir_all(&out_prefix)?;
    } else if let Some(dir) = out_prefix.parent() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir)?;
        }
    }

    let mut writers = vec![];
    for i in 1..=args.n_splits {
        let p = out_prefix.join(format!("{i}.fa"));
        let f = File::create_new(p)?;
        writers.push(BufWriter::new(f))
    }

    let now = std::time::Instant::now();
    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(File::open(args.fa_path)?)?;

    println!("indexed in: {:?}", now.elapsed());

    let mut rng = StdRng::seed_from_u64(args.seed as u64);
    let mut w: Vec<usize> = (0..args.n_splits).collect();

    for i in 1..=index.size {
        let j = i % args.n_splits;
        if j == 0 {
            w.shuffle(&mut rng);
        }
        let seq = index.get(i)?;
        write!(&mut writers[w[j]], "{seq}")?;
    }

    Ok(())
}
