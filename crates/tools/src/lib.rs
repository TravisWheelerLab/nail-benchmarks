use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::Context;

const DIR: &str = env!("CARGO_MANIFEST_DIR");

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
            "missing source data {}; run `make setup` from the repo root",
            path.display()
        );
    }

    Ok(path)
}

pub fn pfam_sto() -> anyhow::Result<PathBuf> {
    data("pfam.sto")
}

// not a make target: built from pfam.sto with hmmbuild
pub fn pfam_hmm() -> anyhow::Result<PathBuf> {
    data("pfam.hmm")
}

pub fn swissprot() -> anyhow::Result<PathBuf> {
    data("swissprot.fa")
}

pub fn mgnify() -> anyhow::Result<PathBuf> {
    data("mgnify")
}
