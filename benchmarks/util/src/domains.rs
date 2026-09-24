//! How a profile's domain hits on one sequence make up its score: which
//! overlap, which run in order, and how much each arrangement adds up to.
//!
//! A port of `classify_direct`, the classifier behind the May 2026 MGnify
//! multi-domain analysis, split into three steps:
//!
//! - [`measure`] scores every arrangement of a pair's domains and takes no
//!   cutoff.
//! - [`classify`] holds a measure against a per-family cutoff and returns the
//!   May analysis's seven categories, exactly as its code did.
//! - [`pattern`] reports which arrangement scores the most, for results that
//!   have no per-family cutoff. It adds [`Pattern::Diffuse`], which has no
//!   counterpart in [`Category`].

use std::collections::{BTreeMap, HashSet};

const OVERLAP_FRAC: f32 = 0.50;
const CORE_COVERAGE: f32 = 0.80;
const EXTEND_RATIO: f32 = 1.20;
const COVERAGE_MID: f32 = 0.50;
const COVERAGE_FULL: f32 = 0.70;
const MIN_FRAG_OVERHANG: u32 = 20;

/// Domains scoring under this take no part in clustering or chaining.
pub const MIN_SCORE: f32 = 4.0;

/// How far under the cutoff a chain or an unordered sum may fall and still
/// count as reaching it.
pub const SLACK: f32 = 4.0;

/// One domain alignment, in model and sequence coordinates, both inclusive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Domain {
    pub score: f32,
    pub hmm_from: u32,
    pub hmm_to: u32,
    pub ali_from: u32,
    pub ali_to: u32,
}

/// The best cluster of two or more overlapping copies of one model region.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cluster {
    pub score: f32,
    pub members: usize,
    /// The model interval covered by at least half the members.
    pub range: (u32, u32),
    /// The length of `range` over the model length.
    pub coverage: f32,
}

/// Everything the classification reads off one pair's domains, before any
/// cutoff is applied.
#[derive(Clone, Debug, PartialEq)]
pub struct Measure {
    /// Domains scoring above zero.
    pub positive: usize,
    /// Positive domains scoring at least [`MIN_SCORE`].
    pub important: usize,
    /// The best single domain's score.
    pub best: f32,
    pub cluster: Option<Cluster>,
    /// The best co-linear chain's score, 0 with no important domain.
    pub chain: f32,
    /// The sum of each cluster's best member, ignoring order.
    pub unordered: f32,
    /// Whether dropping one important domain would leave the rest all
    /// overlapping.
    pub partial: bool,
}

/// How much of the model a cluster's majority range covers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Coverage {
    /// Under half.
    Fragment,
    /// Half to 70%.
    Mid,
    /// 70% or more.
    Full,
}

impl Cluster {
    pub fn band(&self) -> Coverage {
        match self.coverage {
            f if f < COVERAGE_MID => Coverage::Fragment,
            f if f < COVERAGE_FULL => Coverage::Mid,
            _ => Coverage::Full,
        }
    }
}

impl Coverage {
    fn name(self) -> &'static str {
        match self {
            Coverage::Fragment => "fragment",
            Coverage::Mid => "mid",
            Coverage::Full => "full",
        }
    }
}

/// What a pair's score comes to against a per-family cutoff.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Category {
    Strong,
    Clustered(Coverage),
    Ordered,
    Shuffled,
    Weak,
}

impl Category {
    pub const ALL: [Category; 7] = [
        Category::Strong,
        Category::Clustered(Coverage::Fragment),
        Category::Clustered(Coverage::Mid),
        Category::Clustered(Coverage::Full),
        Category::Ordered,
        Category::Shuffled,
        Category::Weak,
    ];

    /// The name the May analysis's tables used.
    pub fn name(self) -> &'static str {
        match self {
            Category::Strong => "strong_match",
            Category::Clustered(Coverage::Fragment) => "clustered:fragment",
            Category::Clustered(Coverage::Mid) => "clustered:mid",
            Category::Clustered(Coverage::Full) => "clustered:full",
            Category::Ordered => "ordered_match",
            Category::Shuffled => "shuffled_match",
            Category::Weak => "weak_match",
        }
    }
}

/// Which arrangement of a pair's domains accounts for the most score, with no
/// cutoff to measure it against.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Pattern {
    /// No domain reaches [`MIN_SCORE`], so the score is spread over hits too
    /// weak to cluster or chain. [`classify`] calls these `weak_match`.
    Diffuse,
    /// No arrangement scores above the best single domain.
    Single,
    /// Overlapping copies of one region of the model.
    Clustered(Coverage),
    /// Domains in the same order in the model and the sequence.
    Ordered,
    /// Domains that only add up with their order ignored.
    Shuffled,
}

impl Pattern {
    pub const ALL: [Pattern; 7] = [
        Pattern::Diffuse,
        Pattern::Single,
        Pattern::Clustered(Coverage::Fragment),
        Pattern::Clustered(Coverage::Mid),
        Pattern::Clustered(Coverage::Full),
        Pattern::Ordered,
        Pattern::Shuffled,
    ];

    pub fn name(self) -> String {
        match self {
            Pattern::Diffuse => "diffuse".to_string(),
            Pattern::Single => "single".to_string(),
            Pattern::Clustered(c) => format!("clustered:{}", c.name()),
            Pattern::Ordered => "ordered".to_string(),
            Pattern::Shuffled => "shuffled".to_string(),
        }
    }
}

/// Measure one pair's domains against a model of `hmm_len` match states, or
/// `None` if no domain scores above zero.
pub fn measure(domains: &[Domain], hmm_len: u32) -> Option<Measure> {
    let mut positive: Vec<Domain> = domains.iter().copied().filter(|d| d.score > 0.0).collect();
    if positive.is_empty() {
        return None;
    }
    positive.sort_by_key(|d| d.hmm_from);

    let best = positive
        .iter()
        .map(|d| d.score)
        .fold(f32::NEG_INFINITY, f32::max);

    let important: Vec<Domain> = positive
        .iter()
        .copied()
        .filter(|d| d.score >= MIN_SCORE)
        .collect();

    let (clusters, doms) = build_clusters(important.clone());
    let score = |c: &[usize]| c.iter().map(|&i| doms[i].score).sum::<f32>();

    let cluster = clusters
        .iter()
        .filter(|c| c.len() >= 2)
        .max_by(|a, b| {
            score(a)
                .partial_cmp(&score(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|c| {
            let (x, y) = majority_range(c, &doms);
            Cluster {
                score: score(c),
                members: c.len(),
                range: (x, y),
                coverage: match y >= x {
                    true => (y - x + 1) as f32 / hmm_len as f32,
                    false => 0.0,
                },
            }
        });

    Some(Measure {
        positive: positive.len(),
        important: important.len(),
        best,
        cluster,
        chain: match important.is_empty() {
            true => 0.0,
            false => chain(&clusters, &doms),
        },
        unordered: clusters
            .iter()
            .map(|c| {
                c.iter()
                    .map(|&i| doms[i].score)
                    .fold(f32::NEG_INFINITY, f32::max)
            })
            .sum(),
        partial: partial_cluster(&important),
    })
}

/// The category a measured pair falls in under a per-family cutoff, and
/// whether it carries the partial-cluster flag.
pub fn classify(m: &Measure, h_cut: f32) -> (Category, bool) {
    if m.best > h_cut {
        return (Category::Strong, false);
    }
    if m.important == 0 {
        return (Category::Weak, false);
    }

    if let Some(c) = m.cluster
        && c.score >= m.chain
    {
        return (Category::Clustered(c.band()), false);
    }

    let category = match () {
        _ if m.chain >= h_cut - SLACK => Category::Ordered,
        _ if m.unordered >= h_cut - SLACK => Category::Shuffled,
        _ => Category::Weak,
    };
    (category, m.partial)
}

/// The arrangement of a measured pair's domains that scores the most, or
/// [`Pattern::Diffuse`] if no domain is strong enough to be arranged.
pub fn pattern(m: &Measure) -> Pattern {
    // the same test classify makes before it looks at any arrangement. with no
    // domain at MIN_SCORE there is no cluster and no chain, and the best
    // single domain would win by default however little of the score it holds
    if m.important == 0 {
        return Pattern::Diffuse;
    }

    // strictly greater, so a tie goes to the simpler arrangement: a chain of
    // one domain is the best single domain, and an unordered sum usually
    // holds the best chain
    let mut best = (Pattern::Single, m.best);
    if let Some(c) = m.cluster
        && c.score > best.1
    {
        best = (Pattern::Clustered(c.band()), c.score);
    }
    if m.chain > best.1 {
        best = (Pattern::Ordered, m.chain);
    }
    if m.unordered > best.1 {
        best = (Pattern::Shuffled, m.unordered);
    }
    best.0
}

// ---

/// The model interval covered by at least half of `members`: the
/// ceil(n/2)-th smallest start to the ceil(n/2)-th largest end.
fn majority_range(members: &[usize], doms: &[Domain]) -> (u32, u32) {
    let n = members.len();
    let k = n.div_ceil(2);
    let mut starts: Vec<u32> = members.iter().map(|&i| doms[i].hmm_from).collect();
    let mut ends: Vec<u32> = members.iter().map(|&i| doms[i].hmm_to).collect();
    starts.sort_unstable();
    ends.sort_unstable();
    (starts[k - 1], ends[n - k])
}

/// Whether two domains share more than half of the shorter one's model span.
fn overlaps(a: &Domain, b: &Domain) -> bool {
    if a.hmm_to < b.hmm_from || b.hmm_to < a.hmm_from {
        return false;
    }
    let span_a = (a.hmm_to - a.hmm_from + 1) as f32;
    let span_b = (b.hmm_to - b.hmm_from + 1) as f32;
    let inter = (a.hmm_to.min(b.hmm_to) - a.hmm_from.max(b.hmm_from) + 1) as f32;
    inter / span_a.min(span_b) > OVERLAP_FRAC
}

/// The connected components of the overlap graph over `h`.
fn components(h: &[usize], doms: &[Domain]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..h.len()).collect();

    fn root(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = root(parent, parent[x]);
        }
        parent[x]
    }

    for i in 0..h.len() {
        for j in (i + 1)..h.len() {
            if overlaps(&doms[h[i]], &doms[h[j]]) {
                let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
                parent[ri.max(rj)] = ri.min(rj);
            }
        }
    }

    // keyed by root rather than hashed, so the order components come out in,
    // and with it every tie below, is the same on every run
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, &d) in h.iter().enumerate() {
        let r = root(&mut parent, i);
        groups.entry(r).or_default().push(d);
    }
    groups.into_values().collect()
}

/// Where in the sequence model position `pos` falls within domain `d`,
/// interpolated linearly.
fn ali_at(pos: u32, d: &Domain) -> u32 {
    let hmm_span = d.hmm_to - d.hmm_from;
    if hmm_span == 0 {
        return d.ali_from;
    }
    // note: a core slice can start before the domain it is cut from, and the
    //       original, built for release, wrapped here rather than panicking.
    //       wrapping explicitly keeps its answers in every build
    let frac = pos.wrapping_sub(d.hmm_from) as f32 / hmm_span as f32;
    d.ali_from
        .wrapping_add((frac * (d.ali_to - d.ali_from) as f32).round() as u32)
}

/// The clusters of `important`, as indices into the returned domains, which
/// hold the originals followed by any pieces split off them.
fn build_clusters(important: Vec<Domain>) -> (Vec<Vec<usize>>, Vec<Domain>) {
    let mut doms = important;
    let mut h: Vec<usize> = (0..doms.len()).collect();
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    while !h.is_empty() {
        let max_score = |c: &[usize]| {
            c.iter()
                .map(|&i| doms[i].score)
                .fold(f32::NEG_INFINITY, f32::max)
        };

        // the largest component, ties to the one with the strongest member;
        // max_by keeps the last of equals, as the original did
        let c = components(&h, &doms)
            .into_iter()
            .max_by(|a, b| {
                a.len().cmp(&b.len()).then_with(|| {
                    max_score(a)
                        .partial_cmp(&max_score(b))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
            })
            .expect("h is not empty");

        let taken: HashSet<usize> = c.iter().copied().collect();
        h.retain(|i| !taken.contains(i));

        if c.len() == 1 {
            clusters.push(c);
            continue;
        }

        let (x, y) = majority_range(&c, &doms);
        let core_len = match y >= x {
            true => (y - x + 1) as f32,
            false => 0.0,
        };

        let (mut core, rest): (Vec<usize>, Vec<usize>) = c.iter().partition(|&&i| {
            let span = (doms[i].hmm_to - doms[i].hmm_from + 1) as f32;
            span > 0.0 && core_len / span > CORE_COVERAGE
        });

        if core.is_empty() {
            clusters.extend(c.into_iter().map(|i| vec![i]));
            continue;
        }

        let max_core = max_score(&core);

        for i in rest {
            if doms[i].score < EXTEND_RATIO * max_core {
                core.push(i);
                continue;
            }

            // much stronger than the core, so only the part over the core
            // joins it, and each long overhang goes back to be clustered
            // again with its share of the score
            let d = doms[i];
            let span = (d.hmm_to - d.hmm_from + 1) as f32;
            let piece = |from: u32, to: u32, ali_from: u32, ali_to: u32| Domain {
                score: d.score * (to - from + 1) as f32 / span,
                hmm_from: from,
                hmm_to: to,
                ali_from,
                ali_to,
            };

            core.push(doms.len());
            doms.push(piece(x, y, ali_at(x, &d), ali_at(y, &d)));

            if x > d.hmm_from && x - d.hmm_from >= MIN_FRAG_OVERHANG {
                h.push(doms.len());
                doms.push(piece(d.hmm_from, x - 1, d.ali_from, ali_at(x - 1, &d)));
            }
            if d.hmm_to > y && d.hmm_to - y >= MIN_FRAG_OVERHANG {
                h.push(doms.len());
                doms.push(piece(y + 1, d.hmm_to, ali_at(y + 1, &d), d.ali_to));
            }
        }

        clusters.push(core);
    }

    (clusters, doms)
}

/// The best total over chains of domains that run in the same order in the
/// model and the sequence, taken across every cluster's members.
fn chain(clusters: &[Vec<usize>], doms: &[Domain]) -> f32 {
    let mut hits: Vec<usize> = clusters.iter().flatten().copied().collect();
    hits.sort_by_key(|&i| doms[i].hmm_from);

    let mut best: Vec<f32> = hits.iter().map(|&i| doms[i].score).collect();
    for j in 1..hits.len() {
        for i in 0..j {
            let (a, b) = (&doms[hits[i]], &doms[hits[j]]);
            let overlap = match a.hmm_to >= b.hmm_from {
                true => (a.hmm_to - b.hmm_from + 1) as f32,
                false => 0.0,
            };
            let min_len =
                ((a.hmm_to - a.hmm_from + 1) as f32).min((b.hmm_to - b.hmm_from + 1) as f32);

            // a few positions of overlap are tolerated as boundary noise
            if overlap <= 10.0 && overlap < 0.1 * min_len && a.ali_to < b.ali_from {
                best[j] = best[j].max(best[i] + b.score);
            }
        }
    }
    best.into_iter().fold(f32::NEG_INFINITY, f32::max)
}

/// Whether every pair of `doms` overlaps.
fn all_overlap(doms: &[Domain]) -> bool {
    (0..doms.len()).all(|i| ((i + 1)..doms.len()).all(|j| overlaps(&doms[i], &doms[j])))
}

/// Whether there are three or more domains and dropping some one of them
/// leaves the rest all overlapping.
fn partial_cluster(doms: &[Domain]) -> bool {
    doms.len() >= 3
        && (0..doms.len()).any(|skip| {
            let rest: Vec<Domain> = doms
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != skip)
                .map(|(_, d)| *d)
                .collect();
            all_overlap(&rest)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    const HMM_LEN: u32 = 100;
    const H_CUT: f32 = 25.0;

    fn d(score: f32, hmm: (u32, u32), ali: (u32, u32)) -> Domain {
        Domain {
            score,
            hmm_from: hmm.0,
            hmm_to: hmm.1,
            ali_from: ali.0,
            ali_to: ali.1,
        }
    }

    fn category(doms: &[Domain], h_cut: f32) -> (Category, bool) {
        classify(&measure(doms, HMM_LEN).expect("a positive domain"), h_cut)
    }

    #[test]
    fn no_positive_domain_is_not_measured() {
        assert_eq!(measure(&[d(0.0, (1, 50), (1, 50))], HMM_LEN), None);
    }

    #[test]
    fn one_domain_over_the_cutoff_is_strong() {
        let doms = [d(30.0, (1, 90), (1, 90)), d(5.0, (1, 20), (200, 220))];
        assert_eq!(category(&doms, H_CUT), (Category::Strong, false));
    }

    #[test]
    fn overlapping_copies_cluster_by_coverage() {
        let full = [d(15.0, (10, 90), (1, 80)), d(14.0, (12, 88), (200, 280))];
        let mid = [d(15.0, (10, 70), (1, 60)), d(14.0, (12, 68), (200, 260))];
        let fragment = [d(12.0, (10, 30), (1, 20)), d(11.0, (12, 28), (200, 220))];

        assert_eq!(
            category(&full, H_CUT),
            (Category::Clustered(Coverage::Full), false)
        );
        assert_eq!(
            category(&mid, H_CUT),
            (Category::Clustered(Coverage::Mid), false)
        );
        assert_eq!(
            category(&fragment, H_CUT),
            (Category::Clustered(Coverage::Fragment), false)
        );
    }

    #[test]
    fn order_decides_ordered_from_shuffled() {
        let ordered = [d(12.0, (1, 40), (1, 40)), d(12.0, (50, 95), (100, 145))];
        let shuffled = [d(12.0, (1, 40), (200, 240)), d(12.0, (50, 95), (10, 55))];

        let m = measure(&ordered, HMM_LEN).unwrap();
        assert_eq!((m.chain, m.unordered), (24.0, 24.0));
        assert_eq!(classify(&m, H_CUT), (Category::Ordered, false));

        let m = measure(&shuffled, HMM_LEN).unwrap();
        assert_eq!((m.chain, m.unordered), (12.0, 24.0));
        assert_eq!(classify(&m, H_CUT), (Category::Shuffled, false));
    }

    #[test]
    fn too_little_is_weak() {
        let scattered = [d(5.0, (1, 40), (1, 40)), d(5.0, (50, 95), (100, 145))];
        let unimportant = [d(3.0, (1, 40), (1, 40))];

        assert_eq!(category(&scattered, H_CUT), (Category::Weak, false));
        assert_eq!(category(&unimportant, H_CUT), (Category::Weak, false));
    }

    #[test]
    fn one_outlier_away_from_a_cluster_is_flagged() {
        // the chain through the outlier outscores the two overlapping copies,
        // so this is ordered; without the outlier the rest would all overlap
        let doms = [
            d(6.0, (10, 50), (10, 50)),
            d(5.0, (12, 48), (300, 340)),
            d(10.0, (80, 95), (100, 115)),
        ];
        assert_eq!(category(&doms, 18.0), (Category::Ordered, true));
    }

    #[test]
    fn pattern_takes_the_arrangement_that_scores_most() {
        let copies = [d(12.0, (1, 40), (1, 40)), d(12.0, (1, 40), (1, 40))];
        let lone = [d(12.0, (1, 40), (1, 40))];
        let ordered = [d(12.0, (1, 40), (1, 40)), d(12.0, (50, 95), (100, 145))];
        let shuffled = [d(12.0, (1, 40), (200, 240)), d(12.0, (50, 95), (10, 55))];
        let full = [d(15.0, (10, 90), (1, 80)), d(14.0, (12, 88), (200, 280))];

        let p = |doms: &[Domain]| pattern(&measure(doms, HMM_LEN).unwrap());
        assert_eq!(p(&lone), Pattern::Single);
        assert_eq!(p(&ordered), Pattern::Ordered);
        assert_eq!(p(&shuffled), Pattern::Shuffled);
        assert_eq!(p(&full), Pattern::Clustered(Coverage::Full));

        // two identical copies cluster to twice one of them
        assert_eq!(p(&copies), Pattern::Clustered(Coverage::Fragment));
    }

    #[test]
    fn many_weak_domains_are_diffuse_whatever_they_add_up_to() {
        let doms: Vec<Domain> = (0..20)
            .map(|i| d(3.5, (10, 40), (1 + 50 * i, 31 + 50 * i)))
            .collect();
        let m = measure(&doms, HMM_LEN).unwrap();

        assert_eq!(pattern(&m), Pattern::Diffuse);
        assert_eq!(classify(&m, H_CUT), (Category::Weak, false));
    }

    #[test]
    fn the_cutoff_changes_the_category_but_not_the_measure() {
        let doms = [d(12.0, (1, 40), (1, 40)), d(12.0, (50, 95), (100, 145))];
        let m = measure(&doms, HMM_LEN).unwrap();

        assert_eq!(classify(&m, 20.0).0, Category::Ordered);
        assert_eq!(classify(&m, 40.0).0, Category::Weak);
        assert_eq!(classify(&m, 10.0).0, Category::Strong);
    }
}
