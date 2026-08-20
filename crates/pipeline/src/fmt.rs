//! Turning numbers into the strings a reader sees.
//!
//! Shared by the table and the progress lines, which disagree about wording but
//! agree about how a duration, a size and a percentage should look.

/// What goes in a cell with no number behind it.
pub(crate) fn dash() -> String {
    "-".to_string()
}

/// Seconds to two places, or a dash.
pub(crate) fn secs(s: Option<f64>) -> String {
    s.map(|s| format!("{s:.2}")).unwrap_or_else(dash)
}

/// `time`'s `%P`: the CPU something burned over the wall clock it took, so a
/// command that kept four cores busy the whole way through reads 400%.
///
/// Truncated rather than rounded, since `time` divides two integers. `time`
/// writes `?%` when there is no clock to divide by; the rest of the table says
/// `-` when it has no number, so that is what goes here.
pub(crate) fn cpu_pct(cpu_s: f64, wall_s: Option<f64>) -> String {
    match wall_s.filter(|wall| *wall > 0.0) {
        Some(wall) => format!("{:.0}%", (cpu_s / wall * 100.0).floor()),
        None => dash(),
    }
}

/// Peak memory in binary units, to about three significant figures so the column
/// reads at a glance: 940KiB, 10.4MiB, 1.02GiB.
pub(crate) fn bytes(kib: i64) -> String {
    if kib < 0 {
        return dash();
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
