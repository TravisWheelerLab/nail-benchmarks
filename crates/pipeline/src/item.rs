//! One thing a step holds, as a sink sees it.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use crate::closure::Closure;
use crate::cmd::Cmd;
use crate::execute::Status;

/// A command or a closure, whichever a step happens to hold.
///
/// A step holds one kind or the other, never both, but a sink almost always
/// wants to treat them the same way: each has a name, a status, fields and
/// tags. The places they genuinely differ are the methods that give back
/// `None` — a closure has no argv to paste, no stderr file to point at, and no
/// exit code, because none of those are things a thread has.
#[derive(Clone, Copy, Debug)]
pub enum Item<'a> {
    Cmd(&'a Cmd),
    Closure(&'a Closure),
}

impl<'a> Item<'a> {
    /// What to call this in a table or on a progress line.
    pub fn label(self) -> String {
        match self {
            Item::Cmd(cmd) => cmd.label(),
            Item::Closure(closure) => closure.label().to_string(),
        }
    }

    pub fn status(self) -> &'a Status {
        match self {
            Item::Cmd(cmd) => cmd.status(),
            Item::Closure(closure) => closure.status(),
        }
    }

    pub fn fields(self) -> &'a BTreeMap<String, String> {
        match self {
            Item::Cmd(cmd) => cmd.fields(),
            Item::Closure(closure) => closure.fields(),
        }
    }

    pub fn tags(self) -> &'a BTreeSet<String> {
        match self {
            Item::Cmd(cmd) => cmd.tags(),
            Item::Closure(closure) => closure.tags(),
        }
    }

    /// The line you could paste into a shell to get the same thing. A closure
    /// has none, run or not.
    pub fn line(self) -> Option<String> {
        match self {
            Item::Cmd(cmd) => Some(cmd.line()),
            Item::Closure(_) => None,
        }
    }

    /// Where its stderr went, for one whose stderr was routed to a file.
    pub fn stderr_path(self) -> Option<&'a Path> {
        match self {
            Item::Cmd(cmd) => cmd.stderr_path(),
            Item::Closure(_) => None,
        }
    }

    /// What it exited with. A closure never exited: it returned, and its status
    /// already says how that went.
    pub fn exit(self) -> Option<i32> {
        match self {
            Item::Cmd(cmd) => cmd.status().timing().map(|t| t.exit),
            Item::Closure(_) => None,
        }
    }
}
