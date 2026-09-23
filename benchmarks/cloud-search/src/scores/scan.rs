//! The domain layer over a results table: which columns hold a hit, and how
//! its score reads.
//!
//! The framing and the field split are `libsail::lines::Rows` and
//! `libsail::tbl::fields`. What is left here is what libsail has no business
//! knowing: hmmer's inclusion column, and the shapes of the target names this
//! repository searches.

use libsail::tbl::HitColumns;
pub use libsail::tbl::fields;
use libsail::tbl::hmmer::HmmerTable;

/// How many fields a line has, which is what a layout is checked against.
pub fn count(line: &[u8]) -> usize {
    line.split(|b| b.is_ascii_whitespace())
        .filter(|field| !field.is_empty())
        .count()
}

/// The query, target and score of one row, at layout `C`'s own indices.
pub fn hit<C: HitColumns>(line: &[u8]) -> Option<[&[u8]; 3]> {
    fields(line, [C::QUERY, C::TARGET, C::SCORE])
}

/// hmmer's `--tblout`, for the one column `libsail`'s layout does not carry:
/// how many of a pair's domains fall inside the inclusion threshold.
///
/// ```text
/// |  0   | 1 |   2  | 3 |  4   |  5  |  6 | ... | 15 | 16 | 17 |   18
/// # target acc query acc e-value score bias ... dom  rep  inc  description
/// ```
const INC: usize = 17;

/// The query, target, score and `inc` of one `--tblout` row.
pub fn hmmer_hit(line: &[u8]) -> Option<[&[u8]; 4]> {
    fields(
        line,
        [
            HmmerTable::QUERY,
            HmmerTable::TARGET,
            HmmerTable::SCORE,
            INC,
        ],
    )
}

/// Whether a row carries the fields layout `C` calls for.
///
/// A table's rows are all one shape, so this is asked of the first row of a
/// file and no others. What it catches is a table of the wrong tool, where
/// reading field 6 as a score would otherwise give a number.
pub fn fits<C: HitColumns>(line: &[u8]) -> bool {
    count(line) >= C::N_COLUMNS
}

/// One score field, or `None` where it is not a number a cutoff can be applied
/// to -- a NaN and an infinity included, since neither compares against one in
/// a way that means anything.
pub fn score(field: &[u8]) -> Option<f32> {
    let n: f32 = std::str::from_utf8(field).ok()?.parse().ok()?;
    n.is_finite().then_some(n)
}

/// The number in an MGnify protein name: `MGYP` and twelve digits.
///
/// Twelve digits is under 2^40, and the numbers order as the names do, so a
/// target turned into one of these sorts a shard without its name being
/// interned, hashed or compared.
pub fn mgyp(name: &[u8]) -> Option<u64> {
    let digits = name.strip_prefix(b"MGYP")?;
    if digits.len() != 12 {
        return None;
    }

    let mut out = 0u64;
    for &byte in digits {
        let digit = byte.wrapping_sub(b'0');
        if digit > 9 {
            return None;
        }

        out = out * 10 + digit as u64;
    }

    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A row of each of the four layouts, as the tool itself wrote it, so the
    /// indices `libsail` declares are held against a real file.
    #[test]
    fn a_hit_is_read_at_its_layout_s_indices() {
        use libsail::tbl::blast::BlastTable;
        use libsail::tbl::hmmer::HmmerDomTable;
        use libsail::tbl::nail::NailTable;

        let nail = b"MGYP005808827855 AAA_30 42     152    22    121   18.1  0.0  9.6e-5 0.066";
        assert!(fits::<NailTable>(nail));
        assert_eq!(
            hit::<NailTable>(nail).unwrap(),
            [&b"AAA_30"[..], b"MGYP005808827855", b"18.1"]
        );

        let mmseqs = b"AAA_30 MGYP005808827855 31.0 100 60 2 3 102 40 139 3.2e-12 45.7";
        assert!(fits::<BlastTable>(mmseqs));
        assert_eq!(
            hit::<BlastTable>(mmseqs).unwrap(),
            [&b"AAA_30"[..], b"MGYP005808827855", b"45.7"]
        );

        // the description field carries spaces, so the row splits into more
        // than the eighteen a layout asks for
        let hmmer = b"MGYP003392859013     -          2_5_RNA_ligase2      PF13563.10     1e-15   54.0   0.2   1.3e-15   53.7   0.2   1.2   1   0   0   1   1   1   1 CR=1 FL=0";
        assert!(fits::<HmmerTable>(hmmer));
        assert_eq!(
            hmmer_hit(hmmer).unwrap(),
            [&b"2_5_RNA_ligase2"[..], b"MGYP003392859013", b"54.0", b"1"]
        );

        let dom = b"MGYP000987338150     -            178 2-oxoacid_dh         PF00198.27   232   9.3e-44  145.5   0.5   1   1   1.1e-46   1.1e-43  145.2   0.5     5   143    42   178    39   178 0.97 FL=0";
        assert!(fits::<HmmerDomTable>(dom));
        assert_eq!(
            hit::<HmmerDomTable>(dom).unwrap(),
            [&b"2-oxoacid_dh"[..], b"MGYP000987338150", b"145.2"]
        );

        // a table of the wrong tool, which is what the field count is for
        assert!(!fits::<HmmerTable>(nail));
        assert!(!fits::<BlastTable>(hmmer_short()));
    }

    /// An hmmer row cut to the fields before the description, which is fewer
    /// than blast's twelve are wide.
    fn hmmer_short() -> &'static [u8] {
        b"MGYP003392859013 - 2_5_RNA_ligase2 PF13563.10 1e-15 54.0"
    }

    #[test]
    fn a_score_is_what_parse_gives() {
        for text in ["98.7", "-3", "0", "1e-30", "-0.0", ".5", "1.", "6.02e23"] {
            assert_eq!(score(text.as_bytes()), text.parse::<f32>().ok());
        }

        for text in ["", "-", "nan", "inf", "-inf", "98.7x", "MGYP1"] {
            assert_eq!(score(text.as_bytes()), None, "{text:?}");
        }
    }

    #[test]
    fn an_mgyp_name_is_its_number() {
        assert_eq!(mgyp(b"MGYP000522683479"), Some(522_683_479));
        assert_eq!(mgyp(b"MGYP000000000000"), Some(0));
        assert_eq!(mgyp(b"MGYP999999999999"), Some(999_999_999_999));
        assert!(mgyp(b"MGYP999999999999").unwrap() < 1 << 40);

        for name in [
            &b"MGYP00052268347"[..],    // eleven digits
            b"MGYP0005226834799",       // thirteen
            b"MGYA000522683479",        // another prefix
            b"MGYP00052268347x",
            b"MGYP",
            b"7tm_1",
            b"",
        ] {
            assert_eq!(mgyp(name), None, "{:?}", String::from_utf8_lossy(name));
        }
    }
}
