use std::{
    fmt::Display,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
};

use anyhow::{anyhow, bail, Context};
use indexmap::IndexMap;

pub const HEADER: &str = "# STOCKHOLM 1.0";

#[derive(Default)]
pub struct StockholmRecord {
    pub id: String,
    pub gf_meta: IndexMap<String, String>,
    pub gs_meta: IndexMap<String, String>,
    pub sequences: IndexMap<String, String>,
}

impl StockholmRecord {
    pub fn parse(lines: Vec<String>) -> anyhow::Result<Self> {
        let mut rec = Self::default();
        let mut lines = lines.iter().filter(|line| !line.is_empty());

        match lines.next() {
            Some(line) => {
                if !line.starts_with("# STOCKHOLM") {
                    bail!("buffer doesn't begin with Stockholm header")
                }
            }
            None => bail!("empty buffer"),
        }

        for line in lines {
            if line.starts_with("#=GF") {
                let mut tokens = line.split_whitespace();
                let key = tokens.nth(1).ok_or(anyhow!("no GF key"))?.to_string();

                if key == "ID" {
                    rec.id = tokens.next().ok_or(anyhow!("no ID value"))?.to_string();
                }

                rec.gf_meta.insert(key, line.clone());
            } else if line.starts_with("#=GS") {
                let mut tokens = line.split_whitespace();
                let name = tokens.nth(1).ok_or(anyhow!("no GS name"))?.to_string();

                rec.gs_meta.insert(name.clone(), line.clone());
                rec.sequences.insert(name, "".to_string());
            } else if line.starts_with("#=GC") || line.starts_with("#=GR") {
                // ignore for now
            } else if !line.starts_with("#") {
                let tokens = line.split_whitespace().collect::<Vec<_>>();
                let name = tokens[0];
                let seq = tokens[1];

                rec.sequences
                    .entry(name.to_string())
                    .and_modify(|s| s.push_str(seq));
            } else {
                bail!("unexpected line format:\n{line}");
            }
        }

        if rec.id.is_empty() {
            bail!("failed to find Stockholm record ID");
        }

        Ok(rec)
    }

    pub fn get(&self, name: &str) -> Option<&String> {
        self.sequences.get(name)
    }
}

impl Display for StockholmRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{HEADER}")?;
        self.gf_meta.values().try_for_each(|s| writeln!(f, "{s}"))?;
        writeln!(f)?;
        self.gs_meta.values().try_for_each(|s| writeln!(f, "{s}"))?;
        writeln!(f)?;
        self.sequences
            .iter()
            .try_for_each(|(n, s)| writeln!(f, "{n} {s}"))?;
        writeln!(f, "//")
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

    pub fn get(&self, name: &str) -> Option<&StockholmRecord> {
        self.records.get(name)
    }

    pub fn get_mut(&mut self, name: &str) -> Option<&mut StockholmRecord> {
        self.records.get_mut(name)
    }

    pub fn write<W: Write>(&self, buf: W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(buf);
        self.records
            .values()
            .try_for_each(|r| writeln!(out, "{r}"))
            .context("Stockholm::write() failed")
    }

    pub fn seq_cnt(&self) -> usize {
        self.records.values().map(|r| r.sequences.len()).sum()
    }
}

/// Copy the records whose `#=GF ID` is in `names` from `src` to `dst`.
///
/// Streams rather than parsing: this is used to cut a subset out of Pfam, which
/// is 500MB, and only the identifier of each block matters.
pub fn subset_by_id(
    src: impl AsRef<Path>,
    names: &std::collections::HashSet<String>,
    dst: impl AsRef<Path>,
) -> anyhow::Result<usize> {
    use std::io::{BufRead, BufReader, BufWriter, Write};

    let src = src.as_ref();
    let mut reader = BufReader::new(
        std::fs::File::open(src)
            .with_context(|| format!("failed to open {}", src.display()))?,
    );
    let mut writer = BufWriter::new(std::fs::File::create(dst.as_ref())?);

    let mut block: Vec<String> = Vec::new();
    let mut id: Option<String> = None;
    let mut kept = 0usize;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix("#=GF ID") {
            id = Some(rest.trim().to_string());
        }

        block.push(line.clone());

        if line.starts_with("//") {
            if id.as_ref().is_some_and(|i| names.contains(i)) {
                for l in &block {
                    writer.write_all(l.as_bytes())?;
                }
                kept += 1;
                if kept == names.len() {
                    break;
                }
            }
            block.clear();
            id = None;
        }
    }

    writer.flush()?;
    Ok(kept)
}

/// Write each record whose `#=GF ID` is in `names` to its own
/// `dst_dir/<id>.sto`, returning how many were written.
///
/// Same single pass as [`subset_by_id`], fanning out instead of concatenating.
/// Useful when a downstream tool takes one alignment at a time.
pub fn explode(
    src: impl AsRef<Path>,
    names: &std::collections::HashSet<String>,
    dst_dir: impl AsRef<Path>,
) -> anyhow::Result<usize> {
    use std::io::{BufRead, BufReader};

    let src = src.as_ref();
    let dst_dir = dst_dir.as_ref();
    std::fs::create_dir_all(dst_dir)?;
    crate::check_file_names(names)?;

    let mut reader = BufReader::new(
        std::fs::File::open(src).with_context(|| format!("failed to open {}", src.display()))?,
    );

    let mut block: Vec<u8> = Vec::new();
    let mut id: Option<String> = None;
    let mut kept = 0usize;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix("#=GF ID") {
            id = Some(rest.trim().to_string());
        }

        block.extend_from_slice(line.as_bytes());

        if line.starts_with("//") {
            if let Some(id) = id.take().filter(|i| names.contains(i)) {
                let path = dst_dir.join(format!("{id}.sto"));
                std::fs::write(&path, &block)
                    .with_context(|| format!("failed to write {}", path.display()))?;
                kept += 1;
                if kept == names.len() {
                    break;
                }
            }
            block.clear();
        }
    }

    Ok(kept)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO: &str = "\
# STOCKHOLM 1.0
#=GF ID alpha
s1 ACGTACGT
s2 ACGTACGT
//
# STOCKHOLM 1.0
#=GF ID beta
s3 WWWWWWWW
//
";

    fn tmp(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("bioio-sto-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.sto");
        std::fs::write(&path, TWO).unwrap();
        path
    }

    #[test]
    fn explode_writes_one_file_per_named_record() {
        let src = tmp("explode");
        let out = src.with_file_name("split");

        let names: std::collections::HashSet<String> = ["beta".to_string()].into_iter().collect();
        let kept = explode(&src, &names, &out).unwrap();

        assert_eq!(kept, 1);
        assert!(!out.join("alpha.sto").exists());

        let text = std::fs::read_to_string(out.join("beta.sto")).unwrap();
        assert!(text.starts_with("# STOCKHOLM 1.0"));
        assert!(text.contains("s3 WWWWWWWW"));
        assert!(text.trim_end().ends_with("//"));
        assert!(!text.contains("alpha"));

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn subset_by_id_keeps_whole_blocks() {
        let src = tmp("subset");
        let dst = src.with_file_name("out.sto");

        let names: std::collections::HashSet<String> = ["alpha".to_string()].into_iter().collect();
        let kept = subset_by_id(&src, &names, &dst).unwrap();

        assert_eq!(kept, 1);
        let text = std::fs::read_to_string(&dst).unwrap();
        assert!(text.contains("s1 ACGTACGT"));
        assert!(!text.contains("beta"));

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }
}
