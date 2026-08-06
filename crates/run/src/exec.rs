use std::fs::{File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

use anyhow::{bail, Context, Result};

use crate::config::Run;

/// Resource usage for one child process, taken from `wait4`.
///
/// This replaces the old `/usr/bin/time -v` wrapper: the numbers come straight
/// out of the kernel rather than being formatted to text and regex-parsed back.
#[derive(Clone, Copy, Debug)]
pub struct Timing {
    pub wall_s: f64,
    pub user_s: f64,
    pub sys_s: f64,
    pub max_rss_kb: i64,
    pub exit: i32,
}

impl Timing {
    /// Aggregate several concurrent children into one logical run: wall clock
    /// is the batch span, CPU time sums, peak RSS is the high-water mark, and
    /// the exit code is the first nonzero one.
    pub fn combine(parts: &[Timing], wall_s: f64) -> Timing {
        Timing {
            wall_s,
            user_s: parts.iter().map(|p| p.user_s).sum(),
            sys_s: parts.iter().map(|p| p.sys_s).sum(),
            max_rss_kb: parts.iter().map(|p| p.max_rss_kb).max().unwrap_or(0),
            exit: parts.iter().map(|p| p.exit).find(|c| *c != 0).unwrap_or(0),
        }
    }
}

/// Where a job's stdout goes. Several tools (lastal, mmseqs convertalis) write
/// their actual results to stdout rather than to a named output file, and
/// psiblast is invoked once per family with its output appended to one table.
#[derive(Clone, Debug)]
pub enum Output {
    Null,
    File(PathBuf),
    Append(PathBuf),
}

#[derive(Clone, Debug)]
pub struct Job {
    pub program: PathBuf,
    pub args: Vec<String>,
    pub stdout: Output,
    pub stderr: Output,
}

impl Job {
    pub fn new(program: impl Into<PathBuf>) -> Self {
        Job {
            program: program.into(),
            args: Vec::new(),
            stdout: Output::Null,
            stderr: Output::Null,
        }
    }

    pub fn arg(mut self, arg: impl Into<String>) -> Self {
        self.args.push(arg.into());
        self
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.args.extend(args.into_iter().map(Into::into));
        self
    }

    pub fn stdout_to(mut self, path: impl Into<PathBuf>) -> Self {
        self.stdout = Output::File(path.into());
        self
    }

    pub fn stdout_append(mut self, path: impl Into<PathBuf>) -> Self {
        self.stdout = Output::Append(path.into());
        self
    }

    /// Capture stderr rather than letting it reach the console. Tools are
    /// chatty (diamond narrates every seed shape), but the output is the only
    /// diagnostic when a run fails, so it is kept until the run succeeds.
    pub fn stderr_to(mut self, path: impl Into<PathBuf>) -> Self {
        self.stderr = Output::Append(path.into());
        self
    }

    /// The full argv as a single line, for the `cmd` column of runs.tsv.
    pub fn display(&self, numa: Option<&Numa>) -> String {
        let mut parts: Vec<String> = Vec::new();
        if let Some(numa) = numa {
            parts.extend(numa.prefix());
        }
        parts.push(self.program.display().to_string());
        parts.extend(self.args.iter().cloned());
        parts.join(" ")
    }

    fn command(&self, numa: Option<&Numa>) -> Result<Command> {
        let mut cmd = match numa {
            Some(numa) => {
                let prefix = numa.prefix();
                let mut cmd = Command::new(&prefix[0]);
                cmd.args(&prefix[1..]);
                cmd.arg(&self.program);
                cmd
            }
            None => Command::new(&self.program),
        };

        cmd.args(&self.args);
        cmd.stdout(stdio(&self.stdout)?);
        cmd.stderr(stdio(&self.stderr)?);

        Ok(cmd)
    }
}

fn stdio(out: &Output) -> Result<Stdio> {
    Ok(match out {
        Output::Null => Stdio::null(),
        Output::File(path) => {
            if let Some(dir) = path.parent() {
                std::fs::create_dir_all(dir)?;
            }
            let file = File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            Stdio::from(file)
        }
        Output::Append(path) => {
            if let Some(dir) = path.parent() {
                std::fs::create_dir_all(dir)?;
            }
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .with_context(|| format!("failed to open {} for append", path.display()))?;
            Stdio::from(file)
        }
    })
}

/// CPU pinning via `numactl`. Only constructed when a config actually asks for
/// a node, so a non-NUMA machine never invokes `numactl` at all.
#[derive(Clone, Debug)]
pub struct Numa {
    node: usize,
    cpus: Vec<usize>,
}

impl Numa {
    /// Pin to the first `threads` CPUs of `node`, mirroring the CPU list that
    /// `set_numa_prefix` in the old util/scripts/run.sh computed via awk.
    pub fn new(node: usize, threads: usize) -> Result<Self> {
        let out = Command::new("numactl")
            .arg("--hardware")
            .output()
            .context("failed to run `numactl --hardware`; is numactl installed?")?;

        if !out.status.success() {
            bail!("`numactl --hardware` failed with status {}", out.status);
        }

        let text = String::from_utf8_lossy(&out.stdout);
        let line = text
            .lines()
            .find(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                f.len() > 3 && f[0] == "node" && f[1] == node.to_string() && f[2] == "cpus:"
            })
            .with_context(|| format!("numa node {node} not present in `numactl --hardware`"))?;

        let cpus: Vec<usize> = line
            .split_whitespace()
            .skip(3)
            .filter_map(|c| c.parse().ok())
            .take(threads)
            .collect();

        if cpus.len() < threads {
            bail!(
                "numa node {node} exposes {} cpus, fewer than the {threads} requested",
                cpus.len()
            );
        }

        Ok(Numa { node, cpus })
    }

    fn prefix(&self) -> Vec<String> {
        let list = self
            .cpus
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",");

        vec![
            "numactl".to_string(),
            format!("--physcpubind={list}"),
            format!("--membind={}", self.node),
        ]
    }
}

/// Run one job to completion, returning its resource usage.
pub fn run(job: &Job, numa: Option<&Numa>) -> Result<Timing> {
    let mut cmd = job.command(numa)?;
    let start = Instant::now();

    let child = cmd
        .spawn()
        .with_context(|| format!("failed to spawn {}", job.program.display()))?;

    let timing = wait_timed(child.id() as libc::pid_t, start)?;
    Ok(timing)
}

/// Run several jobs concurrently and aggregate them into a single Timing.
/// This is the replacement for the `parallel` invocations in run.sh.
pub fn run_concurrent(jobs: &[Job], numa: Option<&Numa>) -> Result<Timing> {
    let start = Instant::now();

    let mut pids = Vec::with_capacity(jobs.len());
    for job in jobs {
        let mut cmd = job.command(numa)?;
        let child = cmd
            .spawn()
            .with_context(|| format!("failed to spawn {}", job.program.display()))?;
        pids.push(child.id() as libc::pid_t);
    }

    let mut parts = Vec::with_capacity(pids.len());
    for pid in pids {
        parts.push(wait_timed(pid, start)?);
    }

    Ok(Timing::combine(&parts, start.elapsed().as_secs_f64()))
}

/// Reap a specific child via wait4 so we get its rusage in isolation.
/// getrusage(RUSAGE_CHILDREN) would accumulate across siblings, which is wrong
/// as soon as we run jobs concurrently.
fn wait_timed(pid: libc::pid_t, start: Instant) -> Result<Timing> {
    let mut status: libc::c_int = 0;
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };

    // std::process::Child has no Drop that reaps, so taking the child with
    // wait4 here does not race with anything in std.
    let rc = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
    if rc < 0 {
        return Err(io::Error::last_os_error()).context("wait4 failed");
    }

    let wall_s = start.elapsed().as_secs_f64();

    let exit = if libc::WIFEXITED(status) {
        libc::WEXITSTATUS(status)
    } else if libc::WIFSIGNALED(status) {
        128 + libc::WTERMSIG(status)
    } else {
        -1
    };

    Ok(Timing {
        wall_s,
        user_s: secs(usage.ru_utime),
        sys_s: secs(usage.ru_stime),
        // ru_maxrss is already in kilobytes on Linux
        max_rss_kb: usage.ru_maxrss,
        exit,
    })
}

fn secs(tv: libc::timeval) -> f64 {
    tv.tv_sec as f64 + tv.tv_usec as f64 / 1_000_000.0
}

/// The single run table per results directory: one row per run, replacing the
/// per-run `.time` and `.summary` files.
///
/// Columns are space-padded to the widest cell in each column and the header
/// is commented out with a dashed separator, matching the layout nail uses for
/// its own `.tbl` output. That means widths are not known until the last row
/// arrives, so rows are kept in memory and the whole file is rewritten after
/// each run. The file on disk therefore stays complete and correctly aligned
/// even if a long run is interrupted partway through.
/// How often the table is rewritten while a run is in progress.
///
/// Column widths are not known until the last row, so every rewrite emits the
/// whole file. Doing that per row is quadratic in bytes written, which a
/// thousand-shard sweep already feels and a per-family calibration pass would
/// make ruinous. Rewriting on a timer instead bounds the loss from a kill to
/// whatever landed in the last few seconds.
const FLUSH_EVERY: std::time::Duration = std::time::Duration::from_secs(5);

pub struct RunsTable {
    path: PathBuf,
    sweep_columns: Vec<String>,
    header: Vec<String>,
    rows: Vec<Vec<String>>,
    last_flush: Instant,
    dirty: bool,
}

impl RunsTable {
    /// Create the table at an explicit path. It does not have to live beside
    /// the hit tables — mgnify keeps it above its results directory so it is
    /// not wiped when a run clears results.
    pub fn create_at(path: impl AsRef<Path>, sweep_columns: Vec<String>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir)?;
        }

        let mut header: Vec<String> = ["name", "tool", "query"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        header.extend(sweep_columns.iter().cloned());
        header.extend(["threads", "target"].iter().map(|s| s.to_string()));
        header.extend(
            ["wall(s)", "user(s)", "sys(s)", "max_rss", "exit", "cmd"]
                .iter()
                .map(|s| s.to_string()),
        );

        let mut table = RunsTable {
            path,
            sweep_columns,
            header,
            rows: Vec::new(),
            last_flush: Instant::now(),
            dirty: false,
        };

        table.flush()?;
        Ok(table)
    }

    pub fn append(
        &mut self,
        run: &Run,
        target: &str,
        timing: &Timing,
        cmd: &str,
    ) -> Result<()> {
        let mut row: Vec<String> = vec![
            run.name.clone(),
            run.tool.clone(),
            run.var_str("query").unwrap_or_else(|| "-".to_string()),
        ];

        // a run only carries the axes its own block swept, so anything from
        // another block's union shows up as "-"
        for col in &self.sweep_columns {
            row.push(
                run.var(col)
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
            );
        }

        row.push(run.threads.to_string());
        row.push(target.to_string());
        row.push(format!("{:.2}", timing.wall_s));
        row.push(format!("{:.2}", timing.user_s));
        row.push(format!("{:.2}", timing.sys_s));
        row.push(format_bytes(timing.max_rss_kb));
        row.push(timing.exit.to_string());
        row.push(cmd.to_string());

        self.rows.push(row);
        self.dirty = true;

        if self.last_flush.elapsed() >= FLUSH_EVERY {
            self.flush()?;
        }

        Ok(())
    }

    /// Write the table out unconditionally. Call this when the run ends, so the
    /// file reflects every row rather than the last timed rewrite.
    pub fn flush(&mut self) -> Result<()> {
        std::fs::write(&self.path, render(&self.header, &self.rows))
            .with_context(|| format!("failed to write {}", self.path.display()))?;

        self.last_flush = Instant::now();
        self.dirty = false;
        Ok(())
    }

    /// Whether rows have been appended since the last write.
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Lay out a commented header, a dashed separator, and the data rows, with
/// every column padded to its widest cell.
fn render(header: &[String], rows: &[Vec<String>]) -> String {
    // the "# " comment marker is absorbed into the first column's width so the
    // header labels stay lined up over the data below them
    let mut head = header.to_vec();
    head[0] = format!("# {}", head[0]);

    let mut widths: Vec<usize> = head.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.chars().count());
        }
    }

    let last = widths.len() - 1;
    let mut sep: Vec<String> = widths.iter().map(|w| "-".repeat(*w)).collect();
    sep[0] = format!("# {}", "-".repeat(widths[0].saturating_sub(2)));
    // cmd is unpadded and hundreds of characters wide; underline the label
    // rather than the whole column
    sep[last] = "-".repeat(header[last].chars().count());

    let mut out = String::new();
    write_row(&mut out, &head, &widths);
    write_row(&mut out, &sep, &widths);
    for row in rows {
        write_row(&mut out, row, &widths);
    }

    out
}

/// Peak memory in binary units, kept to about three significant figures so the
/// column reads at a glance: 940KiB, 10.4MiB, 1.02GiB.
fn format_bytes(kib: i64) -> String {
    if kib < 0 {
        return "-".to_string();
    }

    const STEP: f64 = 1024.0;
    let (value, unit) = match kib as f64 {
        v if v < STEP => (v, "KiB"),
        v if v < STEP * STEP => (v / STEP, "MiB"),
        v => (v / (STEP * STEP), "GiB"),
    };

    if value >= 100.0 {
        format!("{value:.0}{unit}")
    } else if value >= 10.0 {
        format!("{value:.1}{unit}")
    } else {
        format!("{value:.2}{unit}")
    }
}

fn write_row(out: &mut String, cells: &[String], widths: &[usize]) {
    let last = cells.len().saturating_sub(1);
    for (i, cell) in cells.iter().enumerate() {
        if i == last {
            // the final column is cmd, which is long and variable; padding it
            // would only add trailing whitespace
            out.push_str(cell);
        } else {
            let pad = widths[i].saturating_sub(cell.chars().count());
            out.push_str(cell);
            out.extend(std::iter::repeat_n(' ', pad + 1));
        }
    }
    out.push('\n');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timing_captures_cpu_and_exit_code() {
        let job = Job::new("/bin/sh").arg("-c").arg("exit 3");
        let timing = run(&job, None).unwrap();
        assert_eq!(timing.exit, 3);
        assert!(timing.wall_s >= 0.0);
        assert!(timing.max_rss_kb > 0, "expected nonzero peak rss");
    }

    #[test]
    fn signalled_child_reports_128_plus_signal() {
        let job = Job::new("/bin/sh").arg("-c").arg("kill -TERM $$");
        let timing = run(&job, None).unwrap();
        assert_eq!(timing.exit, 128 + libc::SIGTERM);
    }

    #[test]
    fn stdout_can_be_redirected_to_a_file() {
        let dir = std::env::temp_dir().join(format!("bm-exec-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("out.txt");

        let job = Job::new("/bin/sh")
            .arg("-c")
            .arg("echo hello")
            .stdout_to(&out);
        let timing = run(&job, None).unwrap();

        assert_eq!(timing.exit, 0);
        assert_eq!(std::fs::read_to_string(&out).unwrap().trim(), "hello");
        std::fs::remove_dir_all(&dir).ok();
    }

    fn strings(cells: &[&str]) -> Vec<String> {
        cells.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn memory_reads_in_binary_units() {
        assert_eq!(format_bytes(940), "940KiB");
        assert_eq!(format_bytes(10_649), "10.4MiB");
        assert_eq!(format_bytes(107_800), "105MiB");
        assert_eq!(format_bytes(1_066_344), "1.02GiB");
        // exactly on a boundary should step up rather than read as 1024KiB
        assert_eq!(format_bytes(1024), "1.00MiB");
    }

    #[test]
    fn rows_pad_to_the_widest_cell_in_each_column() {
        let header = strings(&["name", "s", "cmd"]);
        let rows = vec![
            strings(&["a", "12.0", "/bin/x --y"]),
            strings(&["longer-name", "5.7", "/bin/z"]),
        ];

        let out = render(&header, &rows);
        let lines: Vec<&str> = out.lines().collect();

        assert_eq!(lines[0], "# name      s    cmd");
        assert_eq!(lines[1], "# --------- ---- ---");
        assert_eq!(lines[2], "a           12.0 /bin/x --y");
        assert_eq!(lines[3], "longer-name 5.7  /bin/z");

        // no line carries trailing padding
        for line in &lines {
            assert_eq!(line.trim_end(), *line, "trailing whitespace in {line:?}");
        }
    }

    #[test]
    fn header_marker_does_not_shift_the_columns() {
        let header = strings(&["name", "tool", "cmd"]);
        let rows = vec![
            strings(&["nail-s12.0.prf", "nail", "/bin/nail"]),
            strings(&["hmmer.seq", "hmmer", "/bin/hmmsearch"]),
        ];

        let out = render(&header, &rows);

        // the "# " marker eats into the first column rather than pushing
        // everything right, so the tool column starts at one offset on every
        // line: header, separator, and data alike
        let tool_starts: Vec<usize> = out
            .lines()
            .map(|line| {
                // skip the marker so header and data lines are measured alike
                let from = if line.starts_with("# ") { 2 } else { 0 };
                let gap = line[from..].find(' ').unwrap() + from;
                line[gap..].find(|c: char| c != ' ').unwrap() + gap
            })
            .collect();

        assert!(
            tool_starts.windows(2).all(|w| w[0] == w[1]),
            "columns not aligned across lines: {tool_starts:?}\n{out}"
        );
    }

    #[test]
    fn concurrent_jobs_sum_cpu_but_not_wall_time() {
        let jobs: Vec<Job> = (0..3)
            .map(|_| {
                Job::new("/bin/sh")
                    .arg("-c")
                    .arg("i=0; while [ $i -lt 200000 ]; do i=$((i+1)); done")
            })
            .collect();

        let timing = run_concurrent(&jobs, None).unwrap();
        assert_eq!(timing.exit, 0);
        // three busy children should burn more cpu than the batch took in wall
        // clock, which is the whole reason we sum user time rather than wall
        assert!(
            timing.user_s + timing.sys_s > timing.wall_s,
            "cpu {} should exceed wall {}",
            timing.user_s + timing.sys_s,
            timing.wall_s
        );
    }
}
