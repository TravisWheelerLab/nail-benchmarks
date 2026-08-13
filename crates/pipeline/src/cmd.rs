use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::execute::Status;
use crate::label;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum Output {
    #[default]
    Null,
    Inherit,
    File(PathBuf),
    Append(PathBuf),
    OnFailure(PathBuf),
}

#[derive(Clone, Debug)]
pub struct Cmd {
    pub(crate) name: Option<String>,
    pub(crate) index: Option<(usize, usize)>,
    pub(crate) argv: Vec<String>,
    pub(crate) env: BTreeMap<String, String>,
    pub(crate) dir: Option<PathBuf>,
    pub(crate) timeout: Option<Duration>,
    pub(crate) stdout: Output,
    pub(crate) stderr: Output,
    pub(crate) fields: BTreeMap<String, String>,
    pub(crate) tags: BTreeSet<String>,
    pub(crate) status: Status,
}

impl Cmd {
    pub fn new(program: impl AsRef<Path>) -> Self {
        Cmd {
            name: None,
            index: None,
            argv: vec![program.as_ref().display().to_string()],
            env: BTreeMap::new(),
            dir: None,
            timeout: None,
            stdout: Output::Null,
            stderr: Output::Null,
            fields: BTreeMap::new(),
            tags: BTreeSet::new(),
            status: Status::NotRun,
        }
    }

    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    pub fn arg(mut self, arg: impl std::fmt::Display) -> Self {
        self.argv.push(arg.to_string());
        self
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: std::fmt::Display,
    {
        self.argv.extend(args.into_iter().map(|a| a.to_string()));
        self
    }

    pub fn path(self, path: impl AsRef<Path>) -> Self {
        self.arg(path.as_ref().display())
    }

    pub fn env(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Self {
        self.env.insert(key.into(), value.to_string());
        self
    }

    pub fn dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.dir = Some(dir.into());
        self
    }

    pub fn timeout(mut self, after: Duration) -> Self {
        self.timeout = Some(after);
        self
    }

    pub fn stdout(mut self, out: Output) -> Self {
        self.stdout = out;
        self
    }

    pub fn stderr(mut self, out: Output) -> Self {
        self.stderr = out;
        self
    }

    pub fn stdout_to(self, path: impl Into<PathBuf>) -> Self {
        self.stdout(Output::File(path.into()))
    }

    pub fn stderr_to(self, path: impl Into<PathBuf>) -> Self {
        self.stderr(Output::File(path.into()))
    }

    pub fn field(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Self {
        self.fields.insert(key.into(), value.to_string());
        self
    }

    pub fn tag(mut self, tag: impl Into<String>) -> Self {
        self.tags.insert(tag.into());
        self
    }

    pub fn label(&self) -> String {
        let index = self.index.map(|(step, cmd)| format!("{step}.{cmd}"));
        label::label(index.as_deref(), self.name.as_deref())
    }

    pub fn status(&self) -> &Status {
        &self.status
    }

    pub fn argv(&self) -> &[String] {
        &self.argv
    }

    pub fn fields(&self) -> &BTreeMap<String, String> {
        &self.fields
    }

    pub fn tags(&self) -> &BTreeSet<String> {
        &self.tags
    }

    pub fn stderr_path(&self) -> Option<&Path> {
        match &self.stderr {
            Output::Null | Output::Inherit => None,
            Output::File(p) | Output::Append(p) | Output::OnFailure(p) => Some(p),
        }
    }

    pub fn line(&self) -> String {
        self.env
            .iter()
            .map(|(key, value)| format!("{key}={}", quote(value)))
            .chain(self.argv.iter().map(|a| quote(a)))
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
