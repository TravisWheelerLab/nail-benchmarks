use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::execute::Status;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum Output {
    #[default]
    Null,
    Inherit,
    File(PathBuf),
    Append(PathBuf),
}

#[derive(Clone, Debug)]
pub struct Cmd {
    pub(crate) name: String,
    pub(crate) argv: Vec<String>,
    pub(crate) stdout: Output,
    pub(crate) stderr: Output,
    pub(crate) labels: BTreeMap<String, String>,
    pub(crate) tags: BTreeSet<String>,
    pub(crate) status: Status,
}

impl Cmd {
    pub fn new(name: impl Into<String>, program: impl AsRef<Path>) -> Self {
        Cmd {
            name: name.into(),
            argv: vec![program.as_ref().display().to_string()],
            stdout: Output::Null,
            stderr: Output::Null,
            labels: BTreeMap::new(),
            tags: BTreeSet::new(),
            status: Status::NotRun,
        }
    }

    pub fn arg(mut self, arg: impl Into<String>) -> Self {
        self.argv.push(arg.into());
        self
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.argv.extend(args.into_iter().map(Into::into));
        self
    }

    pub fn path(self, path: impl AsRef<Path>) -> Self {
        self.arg(path.as_ref().display().to_string())
    }

    pub fn stdout(mut self, out: Output) -> Self {
        self.stdout = out;
        self
    }

    pub fn stderr(mut self, out: Output) -> Self {
        self.stderr = out;
        self
    }

    pub fn label(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Self {
        self.labels.insert(key.into(), value.to_string());
        self
    }

    pub fn tag(mut self, tag: impl Into<String>) -> Self {
        self.tags.insert(tag.into());
        self
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn status(&self) -> &Status {
        &self.status
    }

    pub fn line(&self) -> String {
        self.argv
            .iter()
            .map(|a| quote(a))
            .collect::<Vec<_>>()
            .join(" ")
    }
}

fn quote(arg: &str) -> String {
    const SAFE_PUNCT: &str = "_-./=:+,@%^";

    let plain = !arg.is_empty()
        && arg
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || SAFE_PUNCT.contains(c));

    if plain {
        arg.to_string()
    } else {
        format!("'{}'", arg.replace('\'', r"'\''"))
    }
}
