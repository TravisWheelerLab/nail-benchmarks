//! HMMER3 profile files.
//!
//! A `.hmm` file is a sequence of `//`-terminated blocks, each a header of
//! `KEY value` lines followed by the model itself. Everything here streams:
//! Pfam is 1.6GB, and callers generally want block boundaries or a handful of
//! header fields rather than the models.

use std::collections::HashMap;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context};

/// Per-model statistics needed to turn a bit score into a P-value.
pub struct Gumbel {
    /// Gathering threshold: the curated score above which a hit is a member.
    pub ga_score: f32,
    /// The gathering threshold expressed as a P-value.
    pub ga_p_value: f64,
    /// Location parameter of the forward score distribution.
    pub tau: f64,
    /// Scale parameter of the forward score distribution.
    pub lambda: f64,
}

impl Gumbel {
    /// P-value of a bit score under this model's fitted distribution.
    pub fn p_value(&self, score: f64) -> f64 {
        (-self.lambda * (score - self.tau)).exp()
    }
}

/// Read each model's name, gathering threshold, and forward-score distribution.
pub fn parse_stats(path: impl AsRef<Path>) -> anyhow::Result<HashMap<String, Gumbel>> {
    let path = path.as_ref();
    let reader = BufReader::new(
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );

    let mut names = vec![];
    let mut gathering_thresholds = vec![];
    let mut gumbels = vec![];

    for line in reader.lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix("NAME") {
            names.push(rest.split_whitespace().collect::<String>());
        }

        if let Some(rest) = line.strip_prefix("GA") {
            let x = rest
                .split_whitespace()
                .map(|s| s.parse::<f32>())
                .collect::<anyhow::Result<Vec<_>, _>>()
                .with_context(|| format!("unparseable GA line in {}", path.display()))?;

            gathering_thresholds.push((x[0], x[1]));
        }

        if let Some(rest) = line.strip_prefix("STATS LOCAL FORWARD") {
            let x = rest
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect::<anyhow::Result<Vec<_>, _>>()
                .with_context(|| format!("unparseable STATS line in {}", path.display()))?;

            gumbels.push((x[0], x[1]))
        }
    }

    if names.len() != gathering_thresholds.len() || names.len() != gumbels.len() {
        bail!(
            "{} has {} models but {} GA and {} STATS lines",
            path.display(),
            names.len(),
            gathering_thresholds.len(),
            gumbels.len()
        );
    }

    Ok(names
        .into_iter()
        .enumerate()
        .map(|(i, name)| {
            let ga_score = gathering_thresholds[i].0;
            let tau = gumbels[i].0;
            let lambda = gumbels[i].1;
            let ga_p_value = (-lambda * (ga_score as f64 - tau)).exp();
            (
                name,
                Gumbel {
                    ga_score,
                    ga_p_value,
                    tau,
                    lambda,
                },
            )
        })
        .collect())
}

/// Copy the first `n` complete models of `src` into `dst`, returning their
/// names.
///
/// Completion is counted by `//`, not by `NAME`: stopping at the nth name would
/// truncate that model before its emission lines and terminator.
pub fn subset(src: impl AsRef<Path>, n: usize, dst: impl AsRef<Path>) -> anyhow::Result<HashSet<String>> {
    let src = src.as_ref();
    let mut reader = BufReader::new(
        std::fs::File::open(src).with_context(|| format!("failed to open {}", src.display()))?,
    );
    let mut writer = BufWriter::new(std::fs::File::create(dst.as_ref())?);

    let mut names = HashSet::new();
    let mut complete = 0usize;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix("NAME") {
            names.insert(rest.trim().to_string());
        }

        writer.write_all(line.as_bytes())?;

        if line.starts_with("//") {
            complete += 1;
            if complete == n {
                break;
            }
        }
    }

    writer.flush()?;

    if complete < n {
        bail!(
            "asked for {n} models but {} holds only {complete}",
            src.display()
        );
    }

    Ok(names)
}

/// Write each model whose `NAME` is in `names` to its own `dst_dir/<name>.hmm`,
/// returning how many were written.
///
/// Same single pass as [`subset`], fanning out instead of concatenating, and
/// selecting by name rather than by count. Tools that take one profile per
/// invocation need the models as separate files.
pub fn explode(
    src: impl AsRef<Path>,
    names: &HashSet<String>,
    dst_dir: impl AsRef<Path>,
) -> anyhow::Result<usize> {
    let src = src.as_ref();
    let dst_dir = dst_dir.as_ref();
    std::fs::create_dir_all(dst_dir)?;
    crate::check_file_names(names)?;

    let mut reader = BufReader::new(
        std::fs::File::open(src).with_context(|| format!("failed to open {}", src.display()))?,
    );

    let mut block: Vec<u8> = Vec::new();
    let mut name: Option<String> = None;
    let mut kept = 0usize;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix("NAME") {
            name = Some(rest.trim().to_string());
        }

        block.extend_from_slice(line.as_bytes());

        // blocks are closed by `//`, not by the next NAME: stopping at the name
        // would cut the model off before its emission lines and terminator
        if line.starts_with("//") {
            if let Some(name) = name.take().filter(|n| names.contains(n)) {
                let path = dst_dir.join(format!("{name}.hmm"));
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

    fn tmp(name: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("bioio-hmm-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.hmm");
        std::fs::write(&path, body).unwrap();
        path
    }

    const TWO: &str = "\
HMMER3/f
NAME  alpha
LENG  100
GA    24.60 24.60
STATS LOCAL FORWARD -3.1 0.70
//
HMMER3/f
NAME  beta
LENG  250
GA    21.00 21.00
STATS LOCAL FORWARD -4.2 0.71
//
";

    #[test]
    fn subset_keeps_whole_models() {
        let src = tmp("subset", TWO);
        let dst = src.with_file_name("out.hmm");

        let names = subset(&src, 1, &dst).unwrap();
        let text = std::fs::read_to_string(&dst).unwrap();

        assert_eq!(names.len(), 1);
        assert!(names.contains("alpha"));
        // one complete block: taking the nth NAME would drop the terminator
        assert_eq!(text.matches("//").count(), 1);
        assert!(text.contains("LENG  100"));
        assert!(!text.contains("beta"));

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn subset_rejects_asking_for_too_many() {
        let src = tmp("toomany", TWO);
        let dst = src.with_file_name("out.hmm");

        let err = subset(&src, 5, &dst).unwrap_err().to_string();
        assert!(err.contains("only 2"), "unexpected: {err}");

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn explode_writes_one_file_per_named_model() {
        let src = tmp("explode", TWO);
        let out = src.with_file_name("split");

        let names: HashSet<String> = ["beta".to_string()].into_iter().collect();
        let kept = explode(&src, &names, &out).unwrap();

        assert_eq!(kept, 1);
        assert!(!out.join("alpha.hmm").exists());

        let text = std::fs::read_to_string(out.join("beta.hmm")).unwrap();
        // the whole block travels, terminator included
        assert!(text.starts_with("HMMER3/f"));
        assert!(text.contains("LENG  250"));
        assert!(text.trim_end().ends_with("//"));
        assert!(!text.contains("alpha"));

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn explode_rejects_names_that_escape_the_directory() {
        let src = tmp("escape", TWO);
        let out = src.with_file_name("split");

        let names: HashSet<String> = ["../alpha".to_string()].into_iter().collect();
        let err = explode(&src, &names, &out).unwrap_err().to_string();
        assert!(err.contains("cannot be used as a file name"), "got: {err}");

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn stats_pair_names_with_their_distributions() {
        let src = tmp("stats", TWO);
        let stats = parse_stats(&src).unwrap();

        assert_eq!(stats.len(), 2);
        let alpha = &stats["alpha"];
        assert_eq!(alpha.ga_score, 24.60);
        assert_eq!(alpha.lambda, 0.70);

        // a score at the threshold reproduces the threshold's own p-value
        assert!((alpha.p_value(24.60) - alpha.ga_p_value).abs() < 1e-12);
        // and a higher score is more significant
        assert!(alpha.p_value(40.0) < alpha.ga_p_value);

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }
}
