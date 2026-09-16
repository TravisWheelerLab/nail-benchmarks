//! Where a tool reads and writes, named in a file beside its own source.
//!
//! Nothing in this workspace resolves a location. A tool is told where its
//! inputs are and where its outputs go, and the telling is a `paths.toml` in
//! the tool's own crate directory, one table per label:
//!
//! ```toml
//! [toy]
//! set = "../../store/sets/toy/inputs"
//! run = "../../store/sets/toy/outputs/recall"
//! ```
//!
//! A label is a whole set of paths under one name, so a toy run and a real run
//! differ by a word on the command line. Relative paths mean what they look
//! like they mean: they resolve against the file's own directory, which is why
//! a crate needs no notion of a repository.
//!
//! What may go in one of these files is paths, and for a tool that *makes* a
//! dataset, how much of it to make. Nothing about a search belongs here -- no
//! tools, no flags, no sensitivities, no threads, no name templates. An earlier
//! generation of these files built command lines out of `{}` templates, which
//! made the file the thing that decided what ran and left the real logic in
//! strings nothing type-checked. Every value here is deserialized into a typed
//! field and reaches the code as a `usize` or a `PathBuf`.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use serde::de::DeserializeOwned;

/// What every crate calls its own file.
pub const FILE: &str = "paths.toml";

/// A crate's path file, read.
pub struct File {
    /// What a relative path in the file resolves against.
    dir: PathBuf,
    path: PathBuf,
    labels: toml::Table,
}

impl File {
    /// The file in a crate's own directory. Callers pass
    /// `env!("CARGO_MANIFEST_DIR")`: a crate knowing where its own source lives
    /// is not the same as reaching into another's.
    pub fn open(dir: impl Into<PathBuf>) -> anyhow::Result<File> {
        let dir = dir.into();
        let path = dir.join(FILE);

        let text = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;

        let labels: toml::Table = toml::from_str(&text)
            .with_context(|| format!("failed to parse {}", path.display()))?;

        if labels.is_empty() {
            bail!("{} names no labels", path.display());
        }

        Ok(File { dir, path, labels })
    }

    /// The labels the file names, in the order it names them.
    pub fn labels(&self) -> impl Iterator<Item = &str> {
        self.labels.keys().map(String::as_str)
    }

    /// One label, as whatever shape the caller needs it to be.
    ///
    /// The caller's type decides what a label may contain, so a key that does
    /// not belong is refused here rather than ignored.
    pub fn get<T: DeserializeOwned>(&self, label: &str) -> anyhow::Result<T> {
        let table = self
            .labels
            .get(label)
            .with_context(|| format!("no label {label:?} in {}", self.path.display()))?;

        table.clone().try_into().with_context(|| {
            format!("failed to read label {label:?} from {}", self.path.display())
        })
    }

    /// A path the file named, against the directory the file sits in.
    pub fn at(&self, path: impl AsRef<Path>) -> PathBuf {
        let path = path.as_ref();
        let joined = match path.is_absolute() {
            true => path.to_owned(),
            false => self.dir.join(path),
        };

        tidy(&joined)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// What to print when no label was named: the ones there are, so the next
    /// command can be typed without opening the file.
    pub fn listing(&self, usage: &str) -> String {
        let mut out = format!("labels in {}\n\n", self.path.display());

        for label in self.labels() {
            let _ = writeln!(out, "  {label}");
        }

        let _ = write!(out, "\nusage: {usage}");
        out
    }
}

/// Fold away the `.` and `..` a relative path picked up on the way.
//
// lexical rather than `canonicalize`, which needs the path to exist: most of
// these name something a run is about to create. that makes it wrong for a
// path crossing a symlink, which is a trade worth naming -- these are the
// strings every message prints, and `benchmarks/recall/../../store/sets/toy`
// is not a thing anyone should have to read
fn tidy(path: &Path) -> PathBuf {
    use std::path::Component;

    let mut out = PathBuf::new();

    for part in path.components() {
        match part {
            Component::CurDir => {}
            Component::ParentDir => match out.components().next_back() {
                Some(Component::Normal(_)) => {
                    out.pop();
                }
                _ => out.push(part),
            },
            other => out.push(other),
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;

    #[derive(Deserialize, Debug)]
    #[serde(deny_unknown_fields)]
    struct Paths {
        set: PathBuf,
        run: PathBuf,
    }

    fn file(what: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("util-paths-{}-{what}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(FILE), body).unwrap();
        dir
    }

    const TWO: &str = "\
[toy]
set = \"../sets/toy/inputs\"
run = \"../sets/toy/outputs/recall\"

[full]
set = \"/mnt/scratch/mgy-fixed\"
run = \"../sets/mgy-fixed/outputs/recall\"
";

    #[test]
    fn a_label_becomes_the_callers_type() {
        let dir = file("get", TWO);
        let f = File::open(&dir).unwrap();

        let p: Paths = f.get("toy").unwrap();
        assert_eq!(f.at(&p.set), dir.parent().unwrap().join("sets/toy/inputs"));

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_absolute_path_is_left_alone() {
        let dir = file("abs", TWO);
        let f = File::open(&dir).unwrap();

        let p: Paths = f.get("full").unwrap();
        assert_eq!(f.at(&p.set), PathBuf::from("/mnt/scratch/mgy-fixed"));

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_label_that_is_not_there_names_the_file() {
        let dir = file("unknown", TWO);
        let f = File::open(&dir).unwrap();

        let err = f.get::<Paths>("nope").unwrap_err().to_string();
        assert!(err.contains("nope"), "{err}");
        assert!(err.contains(FILE), "{err}");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_key_that_does_not_belong_is_refused() {
        let dir = file("unknown-key", "[toy]\nset = \"a\"\nrun = \"b\"\nthreads = 8\n");
        let f = File::open(&dir).unwrap();

        let err = f.get::<Paths>("toy").unwrap_err().to_string();
        assert!(err.contains("toy"), "{err}");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn the_dots_a_relative_path_picks_up_are_folded_away() {
        assert_eq!(tidy(Path::new("/a/b/../../store/x")), PathBuf::from("/store/x"));
        assert_eq!(tidy(Path::new("/a/./b")), PathBuf::from("/a/b"));
        // nothing above the root to pop, so it stays as written
        assert_eq!(tidy(Path::new("../x")), PathBuf::from("../x"));
    }

    #[test]
    fn the_listing_names_every_label() {
        let dir = file("listing", TWO);
        let f = File::open(&dir).unwrap();

        let out = f.listing("recall run --in <label>");
        assert!(out.contains("toy"), "{out}");
        assert!(out.contains("full"), "{out}");
        assert!(out.contains("usage: recall run --in <label>"), "{out}");

        std::fs::remove_dir_all(&dir).ok();
    }
}
