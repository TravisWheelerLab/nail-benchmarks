use std::{
    fmt::Display,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use anyhow::Context;
use indexmap::IndexMap;

#[derive(Default, Clone, PartialEq)]
pub struct FastaRecord {
    pub name: String,
    pub extra: String,
    pub sequence: String,
}

impl Display for FastaRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, ">{} {}", self.name, self.extra)?;
        let mut chunks = self
            .sequence
            .as_bytes()
            .chunks(60)
            .map(|c| std::str::from_utf8(c).unwrap());

        if let Some(last) = chunks.next_back() {
            chunks.try_for_each(|c| writeln!(f, "{c}"))?;
            write!(f, "{last}")
        } else {
            Ok(())
        }
    }
}

#[derive(Default)]
pub struct Fasta {
    pub records: IndexMap<String, FastaRecord>,
}

impl Fasta {
    pub fn parse<R: Read>(buf: R) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut fasta = Self::default();
        let mut rec = FastaRecord::default();
        for line in reader.lines() {
            let line = line?;
            if let Some(line) = line.strip_prefix(">") {
                if !rec.name.is_empty() {
                    fasta.records.insert(rec.name.to_string(), rec);
                }
                rec = FastaRecord::default();

                let tokens = line.splitn(2, char::is_whitespace).collect::<Vec<_>>();
                rec.name = tokens[0].to_string();
                rec.extra = tokens[1].to_string();
            } else {
                rec.sequence.push_str(&line)
            }
        }
        fasta.records.insert(rec.name.to_string(), rec);

        Ok(fasta)
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path.as_ref())?;
        Self::parse(f).with_context(|| format!("failed to parse Fasta from: {:?}", path.as_ref()))
    }
}
