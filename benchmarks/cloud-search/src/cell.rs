//! The points the sweep visits.

/// One point on the grid: a pair of pruning thresholds, or the unpruned
/// reference.
///
/// `Full` is nail with `--full-dp`, which skips the cloud stage and fills the
/// whole matrix. It is the ceiling every pruned cell is measured against: the
/// most nail can find off a given seed set, and the longest it can take to find
/// it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell {
    Pruned { a: f32, b: f32 },
    Full,
}

impl Cell {
    /// What this cell's results file is called, and what its step is named.
    pub fn label(self) -> String {
        match self {
            Cell::Pruned { a, b } => format!("A{a:.1}-B{b:.1}"),
            Cell::Full => "full".to_string(),
        }
    }
}

/// Every A against every B, then the unpruned cell on the end.
///
/// Cells where A >= B are in here and are expected to come out identical to
/// each other: A prunes against the best score on the current anti-diagonal and
/// B against the best score anywhere, so a local threshold above the global one
/// never binds. They are left in as a check rather than skipped as waste.
pub fn cells(alphas: &[f32], betas: &[f32]) -> Vec<Cell> {
    let mut out: Vec<Cell> = alphas
        .iter()
        .flat_map(|&a| betas.iter().map(move |&b| Cell::Pruned { a, b }))
        .collect();

    out.push(Cell::Full);
    out
}
