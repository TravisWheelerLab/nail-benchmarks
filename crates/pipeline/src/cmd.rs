use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::execute::Status;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum Output {
    #[default]
    Null,
    Inherit,
    File(PathBuf),
    Append(PathBuf),
    OnFailure(PathBuf),
}

/// Something that can stand in as the value of an option.
///
/// There is no blanket impl over Display, so paths get their own and come out
/// without the quotes a Debug print would add.
pub trait Value {
    fn render(self) -> String;
}

macro_rules! value_via_display {
    ($($t:ty),* $(,)?) => {
        $(impl Value for $t {
            fn render(self) -> String {
                self.to_string()
            }
        })*
    };
}

value_via_display!(
    &str,
    String,
    &String,
    char,
    bool,
    u8,
    u16,
    u32,
    u64,
    usize,
    i8,
    i16,
    i32,
    i64,
    isize,
    f32,
    f64,
    std::fmt::Arguments<'_>,
);

impl Value for &Path {
    fn render(self) -> String {
        self.display().to_string()
    }
}

impl Value for PathBuf {
    fn render(self) -> String {
        self.display().to_string()
    }
}

impl Value for &PathBuf {
    fn render(self) -> String {
        self.display().to_string()
    }
}

/// An option: a flag on its own, or a flag with a value after it.
#[derive(Clone, Debug)]
pub(crate) struct Opt {
    flag: String,
    value: Option<String>,
}

/// A command, built but not run.
///
/// The pieces are kept apart rather than in one argv so that the order they are
/// added in doesn't matter. Several of these tools take their query and target
/// as trailing positionals, and an option tacked on after them would be read as
/// another file.
#[derive(Clone, Debug)]
pub struct Cmd {
    pub(crate) name: Option<String>,
    pub(crate) program: PathBuf,
    /// How many cores this asks for, and the ones it was given. The second is
    /// filled in when the pipeline is built, from the first.
    pub(crate) cores: Option<usize>,
    pub(crate) cpus: Vec<usize>,
    pub(crate) sub: Vec<String>,
    pub(crate) opts: Vec<Opt>,
    pub(crate) positionals: Vec<String>,
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
            program: program.as_ref().to_owned(),
            cores: None,
            cpus: Vec::new(),
            sub: Vec::new(),
            opts: Vec::new(),
            positionals: Vec::new(),
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

    /// Pin this command to `cores` physical cores, whichever ones are going.
    /// Overrides whatever its step asked for.
    pub fn cores(mut self, cores: usize) -> Self {
        self.cores = Some(cores);
        self
    }

    /// A subcommand, like the `search` in `mmseqs search`. Call it more than
    /// once for tools that nest them.
    pub fn sub(mut self, sub: impl Into<String>) -> Self {
        self.sub.push(sub.into());
        self
    }

    /// An option that stands alone, like `--allow-overwrite`.
    pub fn flag(mut self, flag: impl Into<String>) -> Self {
        self.opts.push(Opt {
            flag: flag.into(),
            value: None,
        });
        self
    }

    /// An option and the value that follows it, like `-E 10`.
    pub fn arg(mut self, flag: impl Into<String>, value: impl Value) -> Self {
        self.opts.push(Opt {
            flag: flag.into(),
            value: Some(value.render()),
        });
        self
    }

    /// A positional. These come out in the order they were added, after
    /// everything else.
    pub fn path(mut self, path: impl AsRef<Path>) -> Self {
        self.positionals.push(path.as_ref().display().to_string());
        self
    }

    pub fn env(mut self, key: impl Into<String>, value: impl Value) -> Self {
        self.env.insert(key.into(), value.render());
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

    pub fn field(mut self, key: impl Into<String>, value: impl Value) -> Self {
        self.fields.insert(key.into(), value.render());
        self
    }

    pub fn tag(mut self, tag: impl Into<String>) -> Self {
        self.tags.insert(tag.into());
        self
    }

    /// What to call this command in a table or on the progress line: its name if
    /// it was given one, and otherwise the program it runs, with the path in
    /// front of it dropped. Not unique — the argv column is what tells two
    /// `mkdir`s apart.
    pub fn label(&self) -> String {
        match &self.name {
            Some(name) => name.clone(),
            None => match self.program.file_name() {
                Some(file) => file.to_string_lossy().into_owned(),
                None => self.program.display().to_string(),
            },
        }
    }

    pub fn status(&self) -> &Status {
        &self.status
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

    /// Everything after the program, in the order it gets handed to the shell.
    pub(crate) fn args(&self) -> Vec<String> {
        let mut out = Vec::with_capacity(self.sub.len() + self.opts.len() + self.positionals.len());

        out.extend(self.sub.iter().cloned());

        for opt in &self.opts {
            out.push(opt.flag.clone());
            if let Some(value) = &opt.value {
                out.push(value.clone());
            }
        }

        out.extend(self.positionals.iter().cloned());
        out
    }

    /// What actually gets exec'd, with the pinning wrapper in front of it if
    /// this command was given cpus.
    ///
    /// `program` is left alone rather than being rewritten to `taskset`, so an
    /// unnamed command still labels itself after the thing it runs.
    pub(crate) fn argv(&self) -> (PathBuf, Vec<String>) {
        let mut wrapper = crate::cpu::wrapper(&self.cpus);
        if wrapper.is_empty() {
            return (self.program.clone(), self.args());
        }

        let program = PathBuf::from(wrapper.remove(0));
        let mut args = wrapper;
        args.push(self.program.display().to_string());
        args.extend(self.args());
        (program, args)
    }

    /// The command as a line you can paste into a shell and get the same thing.
    ///
    /// A working directory becomes a subshell so pasting it leaves your own
    /// shell where it was. The redirects sit outside it, because the files are
    /// opened before the child moves anywhere, so a relative one lands in the
    /// same place either way.
    pub fn line(&self) -> String {
        let mut parts: Vec<String> = self
            .env
            .iter()
            .map(|(key, value)| format!("{key}={}", quote(value)))
            .collect();

        let (program, args) = self.argv();
        parts.push(quote(&program.display().to_string()));
        parts.extend(args.iter().map(|a| quote(a)));

        let mut line = parts.join(" ");

        if let Some(dir) = &self.dir {
            line = format!("(cd {} && {line})", quote(&dir.display().to_string()));
        }

        for redirect in [redirect(&self.stdout, ""), redirect(&self.stderr, "2")]
            .into_iter()
            .flatten()
        {
            line.push(' ');
            line.push_str(&redirect);
        }

        line
    }
}

fn redirect(out: &Output, fd: &str) -> Option<String> {
    let (op, path) = match out {
        Output::Inherit => return None,
        Output::Null => (">", "/dev/null".to_string()),
        // OnFailure still writes to the file, it just may not survive the run
        Output::File(p) | Output::OnFailure(p) => (">", quote(&p.display().to_string())),
        Output::Append(p) => (">>", quote(&p.display().to_string())),
    };

    Some(format!("{fd}{op} {path}"))
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
