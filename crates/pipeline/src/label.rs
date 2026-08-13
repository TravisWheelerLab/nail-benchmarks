//! How steps and commands are named in output.
//!
//! One place, so the shape of a name is one decision instead of a format string
//! copied to every site that writes one — and so a caller can take it over later
//! without hunting those sites down.

/// `[1.3](boom)` when both are known, `[1.3]` when there is no name, the bare
/// name when there is no index, and empty when there is neither — which only
/// happens outside a pipeline, since building one hands out the indices.
pub(crate) fn label(index: Option<&str>, name: Option<&str>) -> String {
    match (index, name) {
        (Some(index), Some(name)) => format!("[{index}]({name})"),
        (Some(index), None) => format!("[{index}]"),
        (None, Some(name)) => name.to_string(),
        (None, None) => String::new(),
    }
}

/// The same idea as part of a filename: `3-boom`, or `3` with no name.
pub(crate) fn filename(index: usize, name: Option<&str>) -> String {
    match name {
        Some(name) => format!("{index}-{}", safe(name)),
        None => index.to_string(),
    }
}

/// A name a filesystem will accept. The one place a user's string gets
/// rewritten, because a `/` in a name cannot be part of a path.
fn safe(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || "._-".contains(c) {
                c
            } else {
                '_'
            }
        })
        .collect()
}
