use std::{
    io::{BufRead, BufReader, Read},
    path::Path,
};

use anyhow::{anyhow, Context};
use indexmap::IndexMap;

#[derive(Default)]
pub struct StockholmRecord {
    pub id: String,
    pub sequences: IndexMap<String, String>,
}

impl StockholmRecord {
    pub fn parse(lines: Vec<String>) -> anyhow::Result<Self> {
        let mut rec = Self::default();
        for line in lines {
            if line.is_empty() {
                continue;
            }

            if line.starts_with("#=GF ID") {
                rec.id = line
                    .split_whitespace()
                    .last()
                    .ok_or(anyhow!("no id"))?
                    .to_string();
            } else if line.starts_with("#=GS") {
                let name = line
                    .split_whitespace()
                    .nth(1)
                    .ok_or(anyhow!("no seq name"))?
                    .to_string();

                rec.sequences.insert(name, "".to_string());
            } else if !line.starts_with("#") {
                let tokens = line.split_whitespace().collect::<Vec<_>>();
                let name = tokens[0];
                let seq = tokens[1];

                rec.sequences
                    .entry(name.to_string())
                    .and_modify(|s| s.push_str(seq));
            }
        }
        Ok(rec)
    }
}

#[derive(Default)]
pub struct Stockholm {
    pub records: IndexMap<String, StockholmRecord>,
}

impl Stockholm {
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path.as_ref())
            .with_context(|| format!("failed to open: {}", path.as_ref().to_string_lossy()))?;
        Self::parse(f)
    }

    pub fn parse<R: Read>(buf: R) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut bufs: Vec<Vec<String>> = vec![];
        let mut buf: Vec<String> = vec![];

        for line in reader.lines() {
            // skip malformed lines for now
            let line = line.unwrap_or_default();

            if line.starts_with("//") {
                bufs.push(buf);
                buf = vec![];
            } else {
                buf.push(line);
            }
        }

        let mut sto = Self::default();

        for buf in bufs {
            let rec = StockholmRecord::parse(buf)?;
            sto.records.insert(rec.id.clone(), rec);
        }

        Ok(sto)
    }
}
