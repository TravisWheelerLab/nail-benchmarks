//! Where the programs and the sequence data are.
//!
//! Every benchmark shells out to the same handful of binaries and reads from
//! the same downloads, and neither is looked for on `PATH`: `make` puts the
//! tools in `tools/bin/` and the data in `data/`, and a benchmark reads from
//! here rather than guessing. That is what makes a run reproducible -- whichever
//! `hmmsearch` was built for this repo is the one that ran.
//!
//! A tool accessor checks the binary is there and runs it with `-h` before
//! handing back its path, so a missing or broken install fails by name at the
//! front of a pipeline rather than as a mystery exit code an hour in.
//!
//! [`identity`] is the other half: what a result should record about the binary
//! that produced it. A version string is not enough on its own, since nail can
//! be built from a working tree with `make nail NAIL_SRC=...` and two builds
//! that differ both say `nail 0.7.1`. The hash is of the bytes that actually
//! ran, computed here rather than read out of `tools/installed.tbl`, so it
//! cannot be stale in the one direction that matters.

use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::Context;

const DIR: &str = env!("CARGO_MANIFEST_DIR");

/// The repo root: two directories up from `benchmarks/util`, fixed at compile
/// time. Everything else here is relative to it.
pub fn repo() -> PathBuf {
    PathBuf::from(DIR)
        .parent()
        .and_then(Path::parent)
        .expect("can't find repo root")
        .to_owned()
}

pub fn bin() -> anyhow::Result<PathBuf> {
    Ok(repo().join("tools/bin"))
}

fn tool(name: &str, help: &str) -> anyhow::Result<PathBuf> {
    let path = bin()?.join(name);

    if !path.is_file() {
        anyhow::bail!("no {name} binary at {}", path.display());
    }

    let out = Command::new(&path)
        .arg(help)
        .output()
        .with_context(|| format!("couldn't run {} {help}", path.display()))?;

    if !out.status.success() {
        anyhow::bail!("{} {help} exited {}", path.display(), out.status);
    }

    Ok(path)
}

/// What a tool was, as far as a result needs to record it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Identity {
    pub name: String,
    /// The first version-shaped token the binary prints, or `-`.
    pub version: String,
    /// The first 12 hex digits of the sha256 of the binary, which is what tells
    /// two builds of one version apart. The same digest `make` records.
    pub hash: String,
}

impl std::fmt::Display for Identity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {}", self.name, self.version, self.hash)
    }
}

/// The binary a ledger's `tool` column means.
///
/// The column names the tool a result is attributed to, which for three of
/// them is a suite rather than a program: hmmer is searched with `hmmsearch`,
/// blast with `blastp`, last with `lastal`. The rest are named after their one
/// binary.
pub fn binary_of(tool: &str) -> &str {
    match tool {
        "hmmer" => "hmmsearch",
        "blast" => "blastp",
        "last" => "lastal",
        other => other,
    }
}

/// What the binary at `tools/bin/<name>` is right now.
pub fn identity(tool: &str) -> anyhow::Result<Identity> {
    let name = binary_of(tool);
    let path = bin()?.join(name);
    let bytes = std::fs::read(&path).with_context(|| format!("couldn't read {}", path.display()))?;

    let digest = <sha2::Sha256 as sha2::Digest>::digest(&bytes);
    let hash: String = digest.iter().take(6).map(|b| format!("{b:02x}")).collect();

    Ok(Identity {
        name: name.to_string(),
        version: version_of(&path).unwrap_or_else(|| "-".to_string()),
        hash,
    })
}

/// The first version-shaped token the binary prints, from `--version` or from
/// `-h`. Every tool here answers one of the two, and none answers both the
/// same way, so the shape of the token is what is looked for rather than a
/// per-tool spelling.
fn version_of(path: &Path) -> Option<String> {
    // `version` last and bare: mmseqs answers a subcommand rather than a flag,
    // and a tool that has no such subcommand reads it as a filename, fails,
    // and prints nothing that matches
    let mut fallback = None;

    for flag in ["--version", "-h", "version"] {
        let Ok(out) = Command::new(path).arg(flag).output() else {
            continue;
        };
        if !out.status.success() {
            continue;
        }

        let text = String::from_utf8_lossy(&out.stdout);
        if let Some(found) = text.split_whitespace().find_map(semver) {
            return Some(found);
        }

        // mmseqs answers `version` with a commit and no dots in it, which is
        // its version as much as 3.4 is hmmer's. Kept only if nothing better
        // turns up, since a help text's first word is not a version
        if fallback.is_none()
            && let Some(line) = text.lines().next()
            && let [word] = line.split_whitespace().collect::<Vec<_>>()[..]
            && word.len() <= 64
        {
            fallback = Some(word.to_string());
        }
    }

    fallback
}

/// A token that looks like a version: digits, a dot, then more of either.
fn semver(word: &str) -> Option<String> {
    let trimmed = word.trim_matches(|c: char| !c.is_ascii_alphanumeric());
    let (before, after) = trimmed.split_once('.')?;

    let ok = !before.is_empty()
        && before.bytes().all(|b| b.is_ascii_digit())
        && after.starts_with(|c: char| c.is_ascii_digit());

    ok.then(|| trimmed.to_string())
}

pub fn nail() -> anyhow::Result<PathBuf> {
    tool("nail", "-h")
}

pub fn hmmsearch() -> anyhow::Result<PathBuf> {
    tool("hmmsearch", "-h")
}

pub fn phmmer() -> anyhow::Result<PathBuf> {
    tool("phmmer", "-h")
}

pub fn hmmbuild() -> anyhow::Result<PathBuf> {
    tool("hmmbuild", "-h")
}

pub fn hmmemit() -> anyhow::Result<PathBuf> {
    tool("hmmemit", "-h")
}

pub fn esl_seqstat() -> anyhow::Result<PathBuf> {
    tool("esl-seqstat", "-h")
}

pub fn create_profmark() -> anyhow::Result<PathBuf> {
    tool("create-profmark", "-h")
}

pub fn mmseqs() -> anyhow::Result<PathBuf> {
    tool("mmseqs", "-h")
}

pub fn blastp() -> anyhow::Result<PathBuf> {
    tool("blastp", "-h")
}

pub fn psiblast() -> anyhow::Result<PathBuf> {
    tool("psiblast", "-h")
}

pub fn makeblastdb() -> anyhow::Result<PathBuf> {
    tool("makeblastdb", "-h")
}

pub fn lastal() -> anyhow::Result<PathBuf> {
    tool("lastal", "-h")
}

pub fn lastdb() -> anyhow::Result<PathBuf> {
    tool("lastdb", "-h")
}

pub fn diamond() -> anyhow::Result<PathBuf> {
    tool("diamond", "--help")
}

fn data(name: &str) -> anyhow::Result<PathBuf> {
    let path = repo().join("data").join(name);

    if !path.exists() {
        anyhow::bail!(
            "missing source data {}; run `make data` from the repo root",
            path.display()
        );
    }

    Ok(path)
}

pub fn pfam_sto() -> anyhow::Result<PathBuf> {
    data("pfam.sto")
}

/// Pfam's profiles, built from `pfam.sto`.
///
/// Not routed through `data`: `make data` does not produce it, since building
/// it needs hmmer installed first.
pub fn pfam_hmm() -> anyhow::Result<PathBuf> {
    let path = repo().join("data/pfam.hmm");

    if !path.is_file() {
        anyhow::bail!(
            "missing {}; run `make pfam-hmm` from the repo root",
            path.display()
        );
    }

    Ok(path)
}

pub fn swissprot() -> anyhow::Result<PathBuf> {
    data("swissprot.fa")
}

pub fn mgnify() -> anyhow::Result<PathBuf> {
    data("mgnify")
}

/// The per-family score cutoffs mgy holds every hit against.
///
/// Not routed through `data`: this one is checked in rather than downloaded,
/// so `make data` has nothing to say about it being missing. A calibration run
/// writes a replacement, which is then promoted here by hand.
pub fn mgy_cutoffs() -> anyhow::Result<PathBuf> {
    let path = repo().join("data/mgy-cutoffs.tbl");

    if !path.is_file() {
        anyhow::bail!(
            "missing {}; it is checked in, so a fresh clone should have it",
            path.display()
        );
    }

    Ok(path)
}
