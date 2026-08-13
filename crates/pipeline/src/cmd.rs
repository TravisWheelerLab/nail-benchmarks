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
    /// Write here, then throw the file away unless the command failed and
    /// actually said something. What stderr does by default under a
    /// [`Pipeline`](crate::Pipeline), which picks the path.
    OnFailure(PathBuf),
}

#[derive(Clone, Debug)]
pub struct Cmd {
    pub(crate) name: Option<String>,
    /// Which step this belongs to and where in it, both counting from one.
    /// Handed out by [`PipelineBuilder::build`](crate::PipelineBuilder::build),
    /// so a command that has never been in a pipeline has none.
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

    /// What to call this command in output. Optional: without one it goes by
    /// its number alone.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Anything printable, so a numeric flag needs no `to_string`.
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

    /// A path argument. Paths are not [`Display`](std::fmt::Display), which is
    /// why [`arg`](Self::arg) will not take one.
    pub fn path(self, path: impl AsRef<Path>) -> Self {
        self.arg(path.as_ref().display())
    }

    /// Set a variable for this command, on top of whatever the pipeline itself
    /// was started with.
    pub fn env(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Self {
        self.env.insert(key.into(), value.to_string());
        self
    }

    /// Run from here rather than from wherever the pipeline was started. Worth
    /// setting for anything that drops scratch files in the current directory,
    /// since two of those in one place collide.
    pub fn dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.dir = Some(dir.into());
        self
    }

    /// Give up on this command after `after` and kill it, which comes back as
    /// [`Status::TimedOut`]. Without one it runs as long as it likes.
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

    /// Send stdout to this file, truncating it first. Shorthand for the common
    /// case of [`Output::File`].
    pub fn stdout_to(self, path: impl Into<PathBuf>) -> Self {
        self.stdout(Output::File(path.into()))
    }

    /// Send stderr to this file, truncating it first.
    pub fn stderr_to(self, path: impl Into<PathBuf>) -> Self {
        self.stderr(Output::File(path.into()))
    }

    /// Attach a named value, which shows up as a summary-table column.
    pub fn field(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Self {
        self.fields.insert(key.into(), value.to_string());
        self
    }

    pub fn tag(mut self, tag: impl Into<String>) -> Self {
        self.tags.insert(tag.into());
        self
    }

    /// How this command is named in output: `[1.3](boom)`, or `[1.3]` if it was
    /// never given a name.
    pub fn label(&self) -> String {
        let index = self.index.map(|(step, cmd)| format!("{step}.{cmd}"));
        label::label(index.as_deref(), self.name.as_deref())
    }

    pub fn status(&self) -> &Status {
        &self.status
    }

    /// The program and its arguments, as they are handed to the kernel.
    /// [`line`](Self::line) is this shell-quoted, which is right for pasting and
    /// wrong for anything that wants the pieces.
    pub fn argv(&self) -> &[String] {
        &self.argv
    }

    /// Named values for the summary table, in sorted order.
    pub fn fields(&self) -> &BTreeMap<String, String> {
        &self.fields
    }

    /// Named yes/no facts, in sorted order.
    pub fn tags(&self) -> &BTreeSet<String> {
        &self.tags
    }

    /// Where this command's stderr went, if it went anywhere with a name. Not
    /// called `stderr` because that is already the builder.
    pub fn stderr_path(&self) -> Option<&Path> {
        match &self.stderr {
            Output::Null | Output::Inherit => None,
            Output::File(p) | Output::Append(p) | Output::OnFailure(p) => Some(p),
        }
    }

    /// The command as you would paste it into a shell, environment and all.
    ///
    /// [`dir`](Self::dir) is the one thing missing: there is no prefix for it
    /// the way there is for a variable, so a command that sets one has to be run
    /// from there to match.
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
