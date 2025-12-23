use std::path::PathBuf;

use clap::Parser;
use rayon::prelude::*;

use crate::mmseqs_db::{Db, SplitDb};

#[derive(Parser)]
struct Args {
    #[arg(value_name = "dir/")]
    dir_path: PathBuf,

    #[arg(short, long, value_name = "path")]
    out_path: Option<PathBuf>,

    #[arg(short, value_name = "n", default_value_t = 8)]
    threads: usize,
}

mod mmseqs_db {
    use std::{
        fs::{self, File},
        io::{BufRead, BufReader, Read, Seek},
        path::{Path, PathBuf},
    };

    use anyhow::Context;
    use regex::Regex;

    pub struct Index {
        pub offsets: Vec<u64>,
        pub lengths: Vec<u64>,
    }

    impl Index {
        pub fn from_path<P>(path: P) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
        {
            let reader = BufReader::new(
                File::open(&path)
                    .with_context(|| format!("{}", path.as_ref().to_string_lossy()))?,
            );

            let mut offsets = vec![];
            let mut lengths = vec![];
            for line in reader.lines() {
                let line = line?;

                let tokens = line
                    .split('\t')
                    .map(|s| s.parse::<u64>().expect("failed to parse index line"))
                    .collect::<Vec<_>>();

                offsets.push(tokens[1]);
                // -1 because the index lengths include a null byte
                lengths.push(tokens[2] - 1);
            }

            offsets.shrink_to_fit();
            lengths.shrink_to_fit();

            Ok(Self { offsets, lengths })
        }
    }

    pub struct Db {
        pub index: Index,
        pub path: PathBuf,
    }

    impl Db {
        pub fn from_path<P>(path: P) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
        {
            let path = path.as_ref();
            let index_path = path.with_extension("index");

            Ok(Self {
                index: Index::from_path(index_path)?,
                path: path.to_path_buf(),
            })
        }

        pub fn len(&self) -> usize {
            self.index.offsets.len()
        }

        pub fn get(&self, idx: usize) -> anyhow::Result<String> {
            let offset = self.index.offsets[idx];
            let length = self.index.lengths[idx];

            let mut file = File::open(&self.path)?;
            file.seek(std::io::SeekFrom::Start(offset))?;
            let mut buf = vec![0u8; length as usize];
            file.read_exact(&mut buf)?;

            Ok(String::from_utf8(buf)?)
        }
    }

    pub struct SplitDb {
        pub index: Index,
        pub paths: Vec<PathBuf>,
        pub sizes: Vec<u64>,
    }

    impl SplitDb {
        pub fn from_path<P>(path: P) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
        {
            let path = path.as_ref();
            let index_path = path.with_extension("index");

            // let index_path = dir.join(format!("{name}.index"));
            let dir = path.parent().unwrap();
            let name = path.file_stem().unwrap().to_str().unwrap();

            let re = Regex::new(&format!(r"^{}\.\d+$", regex::escape(name)))?;
            let mut paths = fs::read_dir(dir)?
                .filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    p.file_name()
                        .and_then(|s| s.to_str())
                        .map(|s| re.is_match(s))
                        .unwrap_or(false)
                })
                .collect::<Vec<_>>();

            paths.sort_by_key(|p| {
                p.extension()
                    .and_then(|e| e.to_str())
                    .and_then(|s| s.parse::<usize>().ok())
                    .expect("bad DB split suffix")
            });

            let sizes = paths
                .iter()
                .map(|p| fs::metadata(p).unwrap().len())
                .collect::<Vec<_>>();

            let prefix = sizes
                .into_iter()
                .scan(0u64, |acc, x| {
                    let cur = *acc;
                    *acc += x;
                    Some(cur)
                })
                .collect::<Vec<_>>();

            Ok(Self {
                index: Index::from_path(index_path)?,
                paths,
                sizes: prefix,
            })
        }

        pub fn len(&self) -> usize {
            self.index.offsets.len()
        }

        fn offset_idx(&self, offset: u64) -> usize {
            match self.sizes.binary_search(&offset) {
                Ok(idx) => idx,
                Err(idx) => idx,
            }
        }

        pub fn get(&self, idx: usize) -> anyhow::Result<String> {
            let offset = self.index.offsets[idx];
            let length = self.index.lengths[idx];

            let file_idx = match self.sizes.binary_search(&offset) {
                // this means the binary search found the exact offset
                Ok(idx) => idx,
                // this means the binary search landed between offsets
                Err(idx) => idx.saturating_sub(1),
            };

            let offset = offset - self.sizes[file_idx];

            let mut file = File::open(&self.paths[file_idx])?;
            file.seek(std::io::SeekFrom::Start(offset))?;
            let mut buf = vec![0u8; length as usize];
            file.read_exact(&mut buf)?;

            Ok(String::from_utf8(buf)?)
        }
    }
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let prefilter_db = SplitDb::from_path("./tmp/prefilterDB")?;
    let query_db_header = Db::from_path("./tmp/queryDB_h")?;

    for i in 0..prefilter_db.len() {
        let q = query_db_header.get(i)?;
        let x = prefilter_db.get(i)?;
        println!(
            "{} {}",
            q.trim(),
            x.as_bytes().iter().filter(|&&b| b == b'\n').count()
        );
    }

    // let mut out: Box<dyn Write> = match args.out_path {
    //     Some(p) => Box::new(BufWriter::new(File::create(p)?)),
    //     None => Box::new(BufWriter::new(stdout())),
    // };

    Ok(())
}
