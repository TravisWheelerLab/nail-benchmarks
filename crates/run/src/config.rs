use std::collections::HashSet;
use std::fmt;
use std::path::Path;

use anyhow::{bail, Context, Result};
use indexmap::IndexMap;
use serde::Deserialize;

const DEFAULT_THREADS: usize = 8;

/// A single TOML scalar, as it appears in a sweep axis or a defaults entry.
///
/// Variant order matters: serde tries untagged variants top to bottom, so
/// Int must precede Float to keep `2000` an integer rather than `2000.0`.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum Scalar {
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Scalar::Int(v) => write!(f, "{v}"),
            // sweep values land verbatim in run names (nail-s12.0-ms2000), so a
            // whole-numbered float must keep its decimal point instead of
            // collapsing to "12" the way `{}` would render it
            Scalar::Float(v) if v.fract() == 0.0 => write!(f, "{v:.1}"),
            Scalar::Float(v) => write!(f, "{v}"),
            Scalar::Bool(v) => write!(f, "{v}"),
            Scalar::Str(v) => write!(f, "{v}"),
        }
    }
}

/// A sweep entry: either a fixed value or a list of values to expand over.
#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum Values {
    One(Scalar),
    Many(Vec<Scalar>),
}

impl Values {
    fn as_slice(&self) -> &[Scalar] {
        match self {
            Values::One(s) => std::slice::from_ref(s),
            Values::Many(v) => v,
        }
    }
}

#[derive(Debug, Deserialize)]
pub struct Defaults {
    #[serde(default = "default_threads")]
    pub threads: usize,
    /// Absent means no CPU pinning and no `numactl` invocation at all.
    #[serde(default)]
    pub numa_node: Option<usize>,
    /// Any other key is a free template variable, e.g. `evalue`.
    #[serde(flatten)]
    pub vars: IndexMap<String, Scalar>,
}

impl Default for Defaults {
    fn default() -> Self {
        Self {
            threads: DEFAULT_THREADS,
            numa_node: None,
            vars: IndexMap::new(),
        }
    }
}

fn default_threads() -> usize {
    DEFAULT_THREADS
}

#[derive(Debug, Deserialize)]
pub struct RunBlock {
    pub tool: String,
    pub name: String,
    #[serde(default)]
    pub args: String,
    pub threads: Option<usize>,
    /// When set, the query set is split this many threads per worker and run
    /// as several concurrent processes (HMMER scales poorly past a few threads).
    pub threads_per: Option<usize>,
    /// Every remaining key is a sweep axis.
    #[serde(flatten)]
    pub sweep: IndexMap<String, Values>,
}

#[derive(Debug, Deserialize)]
pub struct Config {
    #[serde(default)]
    pub defaults: Defaults,
    #[serde(rename = "run", default)]
    pub runs: Vec<RunBlock>,
}

/// One concrete invocation, after sweep expansion.
#[derive(Clone, Debug)]
pub struct Run {
    pub name: String,
    pub tool: String,
    pub args: Vec<String>,
    pub threads: usize,
    pub threads_per: Option<usize>,
    /// The concrete sweep assignment for this run, for the runs.tsv columns.
    pub vars: IndexMap<String, Scalar>,
}

impl Run {
    pub fn var(&self, key: &str) -> Option<&Scalar> {
        self.vars.get(key)
    }

    /// A sweep value as a string, for tools that switch on it (e.g. `query`).
    pub fn var_str(&self, key: &str) -> Option<String> {
        self.vars.get(key).map(|v| v.to_string())
    }
}

impl Config {
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read config: {}", path.display()))?;
        toml::from_str(&text).with_context(|| format!("failed to parse config: {}", path.display()))
    }

    /// Load a config whose run blocks live under a name other than `[[run]]`.
    ///
    /// One file can then describe several independent stages that share a
    /// `[defaults]` table — mgnify's calibration has a cheap recruitment sweep
    /// and an exhaustive per-family pass, at different parameterizations.
    pub fn from_path_as(path: impl AsRef<Path>, block: &str) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read config: {}", path.display()))?;

        let mut table: toml::Table = toml::from_str(&text)
            .with_context(|| format!("failed to parse config: {}", path.display()))?;

        let defaults = match table.remove("defaults") {
            Some(value) => value.try_into().with_context(|| {
                format!("failed to parse [defaults] in {}", path.display())
            })?,
            None => Defaults::default(),
        };

        let runs: Vec<RunBlock> = match table.remove(block) {
            Some(value) => value.try_into().with_context(|| {
                format!("failed to parse [[{block}]] in {}", path.display())
            })?,
            None => bail!("{} declares no [[{block}]] blocks", path.display()),
        };

        Ok(Config { defaults, runs })
    }

    /// Every run in the config, with sweeps expanded.
    pub fn expand(&self) -> Result<Vec<Run>> {
        let mut out = Vec::new();
        for block in &self.runs {
            out.extend(
                block
                    .expand(&self.defaults)
                    .with_context(|| format!("in run block {:?}", block.name))?,
            );
        }

        let mut seen = HashSet::new();
        for run in &out {
            if !seen.insert(run.name.clone()) {
                bail!("duplicate run name after expansion: {:?}", run.name);
            }
        }

        Ok(out)
    }

    /// Union of every sweep key across all blocks, minus `query`, which gets a
    /// fixed column of its own because every run has one. These become the
    /// variable middle columns of runs.tsv.
    ///
    /// serde's flatten does not preserve the order keys appear in the TOML, so
    /// these are sorted for a stable, predictable header.
    pub fn sweep_columns(&self) -> Vec<String> {
        let mut cols = Vec::new();
        let mut seen = HashSet::new();
        for block in &self.runs {
            for key in block.sweep.keys() {
                if key != "query" && seen.insert(key.clone()) {
                    cols.push(key.clone());
                }
            }
        }
        cols.sort();
        cols
    }
}

impl RunBlock {
    fn expand(&self, defaults: &Defaults) -> Result<Vec<Run>> {
        let threads = self.threads.unwrap_or(defaults.threads);

        let keys: Vec<&String> = self.sweep.keys().collect();
        let axes: Vec<&[Scalar]> = self.sweep.values().map(Values::as_slice).collect();

        for (key, axis) in keys.iter().zip(&axes) {
            if axis.is_empty() {
                bail!("sweep axis {key:?} is an empty list, which expands to zero runs");
            }
        }

        let mut runs = Vec::new();
        for combo in cartesian(&axes) {
            // defaults first so a sweep axis can shadow a default
            let mut vars = defaults.vars.clone();
            vars.insert("threads".to_string(), Scalar::Int(threads as i64));
            for (key, value) in keys.iter().zip(&combo) {
                vars.insert((*key).clone(), (*value).clone());
            }

            let name = interpolate(&self.name, &vars)?;
            let args = interpolate(&self.args, &vars)?;

            // sweep assignment only; defaults and threads are separate columns
            let assigned = keys
                .iter()
                .zip(&combo)
                .map(|(k, v)| ((*k).clone(), (*v).clone()))
                .collect();

            runs.push(Run {
                name,
                tool: self.tool.clone(),
                args: args.split_whitespace().map(str::to_string).collect(),
                threads,
                threads_per: self.threads_per,
                vars: assigned,
            });
        }

        Ok(runs)
    }
}

fn cartesian<'a>(axes: &[&'a [Scalar]]) -> Vec<Vec<&'a Scalar>> {
    let mut out: Vec<Vec<&Scalar>> = vec![Vec::new()];
    for axis in axes {
        let mut next = Vec::with_capacity(out.len() * axis.len());
        for prefix in &out {
            for value in axis.iter() {
                let mut combo = prefix.clone();
                combo.push(value);
                next.push(combo);
            }
        }
        out = next;
    }
    out
}

/// Substitute `{key}` occurrences. An unknown key is an error rather than a
/// silent passthrough, so a typo in a config surfaces before anything runs.
fn interpolate(template: &str, vars: &IndexMap<String, Scalar>) -> Result<String> {
    let mut out = String::with_capacity(template.len());
    let mut rest = template;

    while let Some(open) = rest.find('{') {
        out.push_str(&rest[..open]);
        let after = &rest[open + 1..];
        let close = after
            .find('}')
            .with_context(|| format!("unclosed '{{' in template {template:?}"))?;

        let key = &after[..close];
        let value = vars.get(key).with_context(|| {
            let known = vars.keys().cloned().collect::<Vec<_>>().join(", ");
            format!("template {template:?} references unknown key {{{key}}}; available: {known}")
        })?;

        out.push_str(&value.to_string());
        rest = &after[close + 1..];
    }

    out.push_str(rest);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(text: &str) -> Config {
        toml::from_str(text).expect("config should parse")
    }

    #[test]
    fn whole_floats_keep_their_decimal_point() {
        assert_eq!(Scalar::Float(12.0).to_string(), "12.0");
        assert_eq!(Scalar::Float(5.7).to_string(), "5.7");
        assert_eq!(Scalar::Int(2000).to_string(), "2000");
    }

    #[test]
    fn sweep_expands_as_a_cartesian_product() {
        let cfg = config(
            r#"
            [defaults]
            threads = 24
            evalue = "1e9"

            [[run]]
            tool = "nail"
            query = ["hmm", "fa"]
            s = [7.5, 12.0]
            max_seqs = [2000]
            name = "nail-s{s}-ms{max_seqs}.{query}"
            args = "--mmseqs-s {s} --mmseqs-max-seqs {max_seqs} -E {evalue}"
            "#,
        );

        let runs = cfg.expand().unwrap();
        assert_eq!(runs.len(), 4);

        let names: Vec<&str> = runs.iter().map(|r| r.name.as_str()).collect();
        assert!(names.contains(&"nail-s12.0-ms2000.hmm"));
        assert!(names.contains(&"nail-s7.5-ms2000.fa"));

        let run = runs.iter().find(|r| r.name.ends_with("s12.0-ms2000.hmm")).unwrap();
        assert_eq!(
            run.args,
            vec!["--mmseqs-s", "12.0", "--mmseqs-max-seqs", "2000", "-E", "1e9"]
        );
        assert_eq!(run.threads, 24);
    }

    #[test]
    fn sweep_columns_are_the_union_in_first_seen_order() {
        let cfg = config(
            r#"
            [[run]]
            tool = "nail"
            query = "hmm"
            s = [12.0]
            name = "nail-{s}"

            [[run]]
            tool = "diamond"
            query = "fa"
            preset = ["fast", "slow"]
            name = "diamond-{preset}"
            "#,
        );

        // query is hoisted to its own fixed column, the rest are sorted
        assert_eq!(cfg.sweep_columns(), vec!["preset", "s"]);
    }

    #[test]
    fn unknown_template_key_is_rejected() {
        let cfg = config(
            r#"
            [[run]]
            tool = "nail"
            s = [12.0]
            name = "nail-{typo}"
            "#,
        );

        let err = cfg.expand().unwrap_err().to_string();
        assert!(err.contains("in run block"), "unexpected error: {err}");
    }

    #[test]
    fn duplicate_run_names_are_rejected() {
        let cfg = config(
            r#"
            [[run]]
            tool = "nail"
            s = [12.0, 12.0]
            name = "nail-{s}"
            "#,
        );

        assert!(cfg.expand().unwrap_err().to_string().contains("duplicate run name"));
    }

    #[test]
    fn numa_is_absent_by_default() {
        let cfg = config(r#"[defaults]
            threads = 24
            "#);
        assert_eq!(cfg.defaults.numa_node, None);
        assert_eq!(cfg.defaults.threads, 24);
    }
}
