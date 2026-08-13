//! Run a list of commands, measure each one, and hand the results to whoever
//! asked for them.
//!
//! A [`Cmd`] is one process: an argv, where its output goes, how it should be
//! described, and — once it has run — a [`Status`] saying how it went. A
//! [`Step`] is a named group of commands that run in order or all at once, and
//! that share a failure policy. A [`Pipeline`] is a list of steps.
//!
//! Output is somebody else's job. A pipeline announces what happened to the
//! [`Sink`]s registered on it and knows nothing more than that: [`Progress`]
//! prints to stderr as commands land, [`Table`] writes the summary table, and
//! anything else is an implementation away.
//!
//! `examples/basic.rs` puts it together.
//!
//! Nothing here knows about shells, config files, or what the commands are for.

mod cmd;
mod execute;
mod label;
mod pipeline;
mod progress;
mod sink;
mod step;
mod table;

pub use cmd::{Cmd, Output};
pub use execute::{Status, Timing, execute};
pub use pipeline::{Pipeline, PipelineBuilder};
pub use progress::Progress;
pub use sink::Sink;
pub use step::{OnError, Step, Strategy};
pub use table::{Headers, Mode, Table};
