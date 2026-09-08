//! Reading the `.time` file a search run outside the harness left behind.
//!
//! A run driven by the harness is timed by `michi`, which calls `wait4` and
//! writes what it measured straight into `manifest.tbl`. A run driven by a
//! shell script on someone else's machine is timed by whatever that machine
//! had, and the six things that could have been are not one format:
//!
//! ```text
//! GNU time -v      Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00
//! POSIX time -p    real 0.01
//! BSD time                 0.01 real         0.00 user         0.00 sys
//! BSD time -l              1212416  maximum resident set size
//! bash builtin     real<TAB>0m0.017s
//! zsh builtin      sleep 1  0.00s user 0.00s system 13% cpu 1.003 total
//! ```
//!
//! Label and value swap sides between GNU and BSD, and the wall clock is
//! written five different ways. So rather than pick the format and then read
//! it, every line is offered to every rule and whatever matches is kept. A
//! label appearing where no rule expects it is ignored rather than an error:
//! the file is someone else's output, and the parts of it this repo reads are
//! a small fraction of what the verbose formats print.
//!
//! Only the wall clock is required. Everything else is what a given format
//! happened to carry.

use std::path::Path;

use anyhow::{Context, bail};

/// What one `.time` file says a command cost.
pub struct Timing {
    pub wall_s: f64,
    pub user_s: Option<f64>,
    pub sys_s: Option<f64>,
    pub cpu_pct: Option<f64>,
    pub max_rss_kb: Option<u64>,

    /// The exit status, where the format reports one. Only GNU's verbose form
    /// does.
    pub exit: Option<i32>,
}

pub fn read(path: &Path) -> anyhow::Result<Timing> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    parse(&text).with_context(|| format!("failed to parse {}", path.display()))
}

/// Every measurement the text carries, in whatever format wrote it.
///
/// A file may hold more than one block: the old harness split a query set,
/// ran the parts under `parallel`, and concatenated their `.time` files. Those
/// blocks overlap in time, so the wall clock is the longest of them rather
/// than their sum, while the processor time they used does add up.
pub fn parse(text: &str) -> anyhow::Result<Timing> {
    let mut found = Found::default();

    for line in text.lines() {
        found.gnu(line);
        found.labelled(line);
        found.bsd(line);
        found.zsh(line);
    }

    let Some(wall_s) = found.wall.iter().copied().fold(None, max_of) else {
        bail!("no wall clock time in it; is this the output of a `time` command?");
    };

    Ok(Timing {
        wall_s,
        user_s: found.user.iter().copied().fold(None, sum_of),
        sys_s: found.sys.iter().copied().fold(None, sum_of),

        // a percentage of one block's wall clock says nothing about several
        // blocks that ran at once
        cpu_pct: match found.cpu.len() {
            1 => Some(found.cpu[0]),
            _ => None,
        },

        max_rss_kb: found.max_rss.iter().copied().max(),

        // any block that failed makes the run a failure, so the first
        // non-zero status is the one that matters
        exit: found.exit.iter().copied().find(|&e| e != 0).or_else(|| {
            // an empty list is "the format never said", not "it succeeded"
            found.exit.first().copied()
        }),
    })
}

fn max_of(acc: Option<f64>, x: f64) -> Option<f64> {
    Some(acc.map_or(x, |a: f64| a.max(x)))
}

fn sum_of(acc: Option<f64>, x: f64) -> Option<f64> {
    Some(acc.unwrap_or(0.0) + x)
}

/// Every value any rule recognized, before they are combined.
#[derive(Default)]
struct Found {
    wall: Vec<f64>,
    user: Vec<f64>,
    sys: Vec<f64>,
    cpu: Vec<f64>,
    max_rss: Vec<u64>,
    exit: Vec<i32>,
}

impl Found {
    /// GNU's verbose form: a label, a colon, then the value.
    fn gnu(&mut self, line: &str) {
        let Some((label, value)) = line.split_once(':') else {
            return;
        };

        let label = label.trim();
        let value = value.trim();

        // the wall clock's own label ends in "(h:mm:ss or m:ss)", so the split
        // above lands inside it and the value is the rest of the line
        if label.starts_with("Elapsed") {
            if let Some((_, rest)) = line.rsplit_once("): ") {
                push(&mut self.wall, duration(rest.trim()));
            }
            return;
        }

        match label {
            "User time (seconds)" => push(&mut self.user, duration(value)),
            "System time (seconds)" => push(&mut self.sys, duration(value)),
            "Percent of CPU this job got" => {
                push(&mut self.cpu, value.trim_end_matches('%').parse().ok())
            }
            "Maximum resident set size (kbytes)" => {
                push(&mut self.max_rss, value.parse().ok())
            }
            "Exit status" => push(&mut self.exit, value.parse().ok()),
            _ => {}
        }
    }

    /// The POSIX form and bash's builtin: the label first, then one value.
    ///
    /// `real 0.01` and `real\t0m0.017s` differ only in how the duration is
    /// written, which [`duration`] settles.
    fn labelled(&mut self, line: &str) {
        let mut fields = line.split_whitespace();
        let (Some(label), Some(value), None) = (fields.next(), fields.next(), fields.next())
        else {
            return;
        };

        match label {
            "real" => push(&mut self.wall, duration(value)),
            "user" => push(&mut self.user, duration(value)),
            "sys" => push(&mut self.sys, duration(value)),
            _ => {}
        }
    }

    /// BSD's forms: the value first, then the label.
    ///
    /// The default form puts all three pairs on one line; `-l` adds a block of
    /// `getrusage` lines whose labels are several words long.
    fn bsd(&mut self, line: &str) {
        // zsh writes value-then-label pairs too, so without this both rules
        // match its line and every value on it is counted twice
        if is_zsh(line) {
            return;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();

        // maximum resident set size, whose label runs to the end of the line.
        // macOS reports it in bytes where every other producer uses kilobytes
        if let [value, rest @ ..] = &fields[..]
            && rest.join(" ") == "maximum resident set size"
        {
            push(&mut self.max_rss, value.parse::<u64>().ok().map(|b| b / 1024));
            return;
        }

        for pair in fields.chunks(2) {
            let [value, label] = pair else { continue };

            match *label {
                "real" => push(&mut self.wall, duration(value)),
                "user" => push(&mut self.user, duration(value)),
                "sys" => push(&mut self.sys, duration(value)),
                _ => {}
            }
        }
    }

    /// zsh's builtin, which puts the whole thing on one line and calls the
    /// wall clock `total`.
    fn zsh(&mut self, line: &str) {
        if !is_zsh(line) {
            return;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();

        for pair in fields.windows(2) {
            let [value, label] = pair else { continue };

            match *label {
                "total" => push(&mut self.wall, duration(value)),
                "user" => push(&mut self.user, duration(value)),
                "system" => push(&mut self.sys, duration(value)),
                "cpu" => push(&mut self.cpu, value.trim_end_matches('%').parse().ok()),
                _ => {}
            }
        }
    }
}

/// zsh's builtin is the only producer that names the wall clock `total`, which
/// is enough to tell its line from BSD's.
fn is_zsh(line: &str) -> bool {
    line.split_whitespace().last() == Some("total")
}

fn push<T>(out: &mut Vec<T>, value: Option<T>) {
    if let Some(value) = value {
        out.push(value);
    }
}

/// One duration, in any of the ways the six formats write it.
///
/// ```text
/// 12.34      seconds
/// 0.00s      seconds, zsh's suffix
/// 1:02.5     minutes and seconds
/// 1:02:03    hours, minutes and seconds
/// 0m0.017s   bash's composite
/// ```
fn duration(s: &str) -> Option<f64> {
    let s = s.trim();

    // bash's builtin, which is the only form using a unit letter as a
    // separator rather than a suffix
    if let Some((minutes, seconds)) = s.split_once('m')
        && let Some(seconds) = seconds.strip_suffix('s')
    {
        let minutes: f64 = minutes.parse().ok()?;
        let seconds: f64 = seconds.parse().ok()?;
        return Some(minutes * 60.0 + seconds);
    }

    // h:mm:ss and m:ss both count upwards in 60s from the right, so the
    // number of parts does not have to be known
    let mut total = 0.0;
    for (i, part) in s.trim_end_matches('s').rsplit(':').enumerate() {
        let part: f64 = part.parse().ok()?;
        total += part * 60_f64.powi(i as i32);
    }

    Some(total)
}
