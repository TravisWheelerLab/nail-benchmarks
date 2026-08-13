//! Run a list of commands, measure each one, and write down what happened.
//!
//! The whole model is three types. A [`Cmd`] is one process: an argv, where its
//! output goes, how it should be described in the table, and — once it has run —
//! what it cost. A [`Step`] is a named batch or sequence of commands sharing a
//! failure policy. A [`Pipeline`] is a list of steps plus somewhere to put the
//! [`Report`], which is those same steps handed back with their timings filled
//! in.
//!
//! `examples/basic.rs` puts the three together.
//!
//! Steps always run one after another. Within a step, commands run in order or
//! all at once, depending on how the step was built. Nothing here knows about
//! shells, config files, or what the commands are for.

mod cmd;
mod execute;
mod pipeline;
mod step;
mod table;

pub use cmd::{Cmd, Output};
pub use execute::{Status, Timing, execute};
pub use pipeline::Pipeline;
pub use step::{OnError, Step};
pub use table::Report;
