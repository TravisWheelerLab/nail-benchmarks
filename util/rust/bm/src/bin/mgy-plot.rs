use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use clap::{Parser, Subcommand};

// | GA |    GA     |  nail  |    nail   |  nail  |   nail    | mmseqs |  mmseqs   | hmmer |  hmmer    | dom | dom | sig |
// | sc |  P-value  | cld sc |cld P-value|   sc   |  P-value  |   sc   |  P-value  |  sc   |  P-value  | sum | max | dom | dom scores
//  ---- ----------- -------- ----------- -------- ----------- -------- ----------- ------- ----------- ----- ----- ----- ------------

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    Ga {
        #[arg(value_name = "hits.tbl")]
        tbl_path: PathBuf,
    },
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().cmd {
        Cmd::Ga { tbl_path } => {
            max_seqs_info(tbl_path)?;
        }
    }

    Ok(())
}

#[derive(Default)]
struct GaHitsRecord {
    hmmer_hits: usize,
    nail_hits: usize,
}

fn max_seqs_info(path: impl AsRef<Path>) -> anyhow::Result<()> {
    let reader = BufReader::new(File::open(path)?);

    let mut map: HashMap<String, GaHitsRecord> = HashMap::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            continue;
        }

        let tokens = line.split_whitespace().collect::<Vec<_>>();

        let query = tokens[0];
        let ga_sc = tokens[2].parse::<f32>()?;
        let nail_sc = tokens[6].parse::<f32>().ok();
        let hmmer_sc = tokens[10].parse::<f32>().ok();

        let record = match map.get_mut(query) {
            Some(r) => r,
            None => {
                map.insert(query.to_string(), GaHitsRecord::default());
                map.get_mut(query).unwrap()
            }
        };

        if let Some(sc) = hmmer_sc {
            record.hmmer_hits += (sc >= ga_sc) as usize
        }

        if let Some(sc) = nail_sc {
            record.nail_hits += (sc >= ga_sc) as usize
        }
    }

    let mut hmmer_tot = 0;
    let mut nail_tot = 0;

    for record in map.values() {
        hmmer_tot += record.hmmer_hits;
        nail_tot += record.nail_hits;
    }

    println!("{}", nail_tot as f32 / hmmer_tot as f32);

    Ok(())
}
