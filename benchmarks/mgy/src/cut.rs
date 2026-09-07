//! Cutting Pfam down to a query set, as both an hmm file and its alignments.
//!
//! Two shapes of cut, each in both formats: `subset_*` takes a prefix into one
//! file, and `explode_*` fans records out one file per family. Only this
//! benchmark cuts Pfam up, so they live here rather than in `bench`.

use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, bail};
use libsail::collection::{Indexable, Iterable};
use libsail::index::build_index;
use libsail::parse::Parse;
use libsail::seq::p7hmm::IndexedHmm;
use libsail::seq::stockholm::{StockholmDelimiter, StockholmParser, StockholmRecord};
use libsail::source::{FileSource, Source};

/// Copy the first `n` models of `src` into `dst`, returning their names.
pub fn subset_hmm(
    src: impl AsRef<Path>,
    n: usize,
    dst: impl AsRef<Path>,
) -> anyhow::Result<HashSet<String>> {
    let src = src.as_ref();
    let models =
        IndexedHmm::open(src).with_context(|| format!("failed to index {}", src.display()))?;

    if models.len() < n {
        bail!(
            "asked for {n} models but {} holds only {}",
            src.display(),
            models.len()
        );
    }

    let dst = dst.as_ref();
    let mut writer = BufWriter::new(
        std::fs::File::create(dst)
            .with_context(|| format!("failed to create {}", dst.display()))?,
    );

    let mut names = HashSet::with_capacity(n);

    // take before iter, so only the models that travel are parsed
    for model in models.take(n).iter() {
        names.insert(model.header.name.clone());
        model.write_to(&mut writer)?;
    }

    writer.flush()?;
    Ok(names)
}

/// Write each model whose `NAME` is in `names` to
/// `dst_dir/<name>/query.hmm`, returning how many were written.
///
/// One directory per family holding a `query.hmm` and a `query.sto`, which is
/// the shape a ladder rung's query directory has and what the tools that take
/// one family at a time want.
pub fn explode_hmm(
    src: impl AsRef<Path>,
    names: &HashSet<String>,
    dst_dir: impl AsRef<Path>,
) -> anyhow::Result<usize> {
    let src = src.as_ref();
    let dst_dir = dst_dir.as_ref();
    std::fs::create_dir_all(dst_dir)?;
    check_file_names(names)?;

    let models =
        IndexedHmm::open(src).with_context(|| format!("failed to index {}", src.display()))?;

    let mut kept = 0usize;

    for model in models.iter() {
        if !names.contains(&model.header.name) {
            continue;
        }

        write_one(&dst_dir.join(&model.header.name).join("query.hmm"), |w| {
            model.write_to(w).map_err(Into::into)
        })?;

        kept += 1;
        if kept == names.len() {
            break;
        }
    }

    Ok(kept)
}

/// Copy the alignments whose `#=GF ID` is in `names` from `src` to `dst`,
/// returning how many were written.
pub fn subset_sto(
    src: impl AsRef<Path>,
    names: &HashSet<String>,
    dst: impl AsRef<Path>,
) -> anyhow::Result<usize> {
    let dst = dst.as_ref();
    let mut writer = BufWriter::new(
        std::fs::File::create(dst)
            .with_context(|| format!("failed to create {}", dst.display()))?,
    );

    let kept = for_each_named(src.as_ref(), names, |rec| {
        rec.write_to(&mut writer)?;
        Ok(())
    })?;

    writer.flush()?;
    Ok(kept)
}

/// Write each alignment whose `#=GF ID` is in `names` to
/// `dst_dir/<id>/query.sto`, returning how many were written.
///
/// The same per-family directory [`explode_hmm`] writes into.
pub fn explode_sto(
    src: impl AsRef<Path>,
    names: &HashSet<String>,
    dst_dir: impl AsRef<Path>,
) -> anyhow::Result<usize> {
    let dst_dir = dst_dir.as_ref();
    std::fs::create_dir_all(dst_dir)?;
    check_file_names(names)?;

    for_each_named(src.as_ref(), names, |rec| {
        let id = rec.id().expect("named by the walk that selected it");

        write_one(&dst_dir.join(id).join("query.sto"), |w| {
            rec.write_to(w).map_err(Into::into)
        })
    })
}

/// Hand `f` each alignment of `src` whose `#=GF ID` is in `names`, stopping
/// once every name has been seen.
//
// framed and read through libsail, but parsed a range at a
// time rather than through IndexedStockholm: a record that
// will not parse has to be repaired and parsed again, and
// Indexable::get drops both the bytes and the reason -- it
// answers None, which every combinator reads as a broken
// contract and panics on
fn for_each_named<F>(src: &Path, names: &HashSet<String>, mut f: F) -> anyhow::Result<usize>
where
    F: FnMut(&StockholmRecord) -> anyhow::Result<()>,
{
    let source =
        FileSource::open(src).with_context(|| format!("failed to open {}", src.display()))?;
    let offsets = build_index::<_, StockholmDelimiter>(source.file())
        .with_context(|| format!("failed to index {}", src.display()))?;

    let mut kept = 0usize;

    for (n, offset) in offsets.iter().enumerate() {
        let bytes = source.range(offset.start, offset.n_bytes)?;

        let rec = match StockholmParser::parse(&bytes) {
            Ok(rec) => rec,
            Err(_) => {
                let mut repaired = bytes.into_owned();
                util::repair_utf8(&mut repaired);

                StockholmParser::parse(&repaired)
                    .with_context(|| format!("record {n} of {} will not parse", src.display()))?
            }
        };

        if !rec.id().is_some_and(|id| names.contains(id)) {
            continue;
        }

        f(&rec)?;

        kept += 1;
        if kept == names.len() {
            break;
        }
    }

    Ok(kept)
}

/// Write one record to its own file, making the directory it goes in.
fn write_one<F>(path: &Path, f: F) -> anyhow::Result<()>
where
    F: FnOnce(&mut BufWriter<std::fs::File>) -> anyhow::Result<()>,
{
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("failed to create {}", dir.display()))?;
    }

    let mut writer = BufWriter::new(
        std::fs::File::create(path)
            .with_context(|| format!("failed to write {}", path.display()))?,
    );

    f(&mut writer)?;
    writer.flush()?;
    Ok(())
}

/// Reject record identifiers that cannot safely become file names.
///
/// The `explode` functions turn identifiers into paths, and identifiers come
/// out of whatever file was handed in. A separator or `..` would put the output
/// somewhere the caller did not ask for.
fn check_file_names(names: &HashSet<String>) -> anyhow::Result<()> {
    for name in names {
        if name.is_empty()
            || name == ".."
            || name.contains(std::path::MAIN_SEPARATOR)
            || name.contains('/')
        {
            bail!("record identifier {name:?} cannot be used as a file name");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use util::profile;

    fn tmp(name: &str, body: &[u8]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("mgy-cut-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, body).unwrap();
        path
    }

    fn two_models() -> String {
        format!("{}{}", profile("alpha", 4), profile("beta", 6))
    }

    const TWO_ALIGNMENTS: &str = "\
# STOCKHOLM 1.0
#=GF ID alpha
s1 ACGTACGT
s2 ACGTACGT
//
# STOCKHOLM 1.0
#=GF ID beta
s3 WWWWWWWW
//
";

    fn ids_and_rows(path: &Path) -> Vec<(String, Vec<String>)> {
        let source = FileSource::open(path).unwrap();
        let offsets = build_index::<_, StockholmDelimiter>(source.file()).unwrap();

        offsets
            .iter()
            .map(|o| {
                let bytes = source.range(o.start, o.n_bytes).unwrap();
                let rec = StockholmParser::parse(&bytes).unwrap();

                let rows = rec
                    .names
                    .iter()
                    .zip(&rec.seqs)
                    .map(|(name, seq)| format!("{name} {}", String::from_utf8_lossy(seq)))
                    .collect();

                (rec.id().unwrap().to_string(), rows)
            })
            .collect()
    }

    #[test]
    fn subset_hmm_keeps_whole_models() {
        let src = tmp("subset.hmm", two_models().as_bytes());
        let dst = src.with_file_name("out.hmm");

        let names = subset_hmm(&src, 1, &dst).unwrap();

        assert_eq!(names.len(), 1);
        assert!(names.contains("alpha"));

        // one whole model, readable as one: a truncated block would not index
        let back = IndexedHmm::open(&dst).unwrap();
        assert_eq!(back.len(), 1);

        let rec = back.cloned(0).unwrap();
        assert_eq!(rec.header.name, "alpha");
        assert_eq!(rec.header.leng, 4);

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn subset_hmm_rejects_asking_for_too_many() {
        let src = tmp("toomany.hmm", two_models().as_bytes());
        let dst = src.with_file_name("out.hmm");

        let err = subset_hmm(&src, 5, &dst).unwrap_err().to_string();
        assert!(err.contains("only 2"), "unexpected: {err}");

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn explode_hmm_writes_one_directory_per_named_model() {
        let src = tmp("explode.hmm", two_models().as_bytes());
        let out = src.with_file_name("split");

        let names: HashSet<String> = ["beta".to_string()].into_iter().collect();
        assert_eq!(explode_hmm(&src, &names, &out).unwrap(), 1);
        assert!(!out.join("alpha").exists());

        let back = IndexedHmm::open(out.join("beta").join("query.hmm")).unwrap();
        assert_eq!(back.len(), 1);

        let rec = back.cloned(0).unwrap();
        assert_eq!(rec.header.name, "beta");
        assert_eq!(rec.header.leng, 6);

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    /// Both formats land in the same per-family directory, which is what makes
    /// it a query set rather than two parallel piles.
    #[test]
    fn both_explodes_fill_one_directory_per_family() {
        let hmm = tmp("pair.hmm", two_models().as_bytes());
        let sto = hmm.with_file_name("pair.sto");
        std::fs::write(&sto, TWO_ALIGNMENTS).unwrap();

        let out = hmm.with_file_name("queries");
        let names: HashSet<String> = ["alpha".to_string(), "beta".to_string()]
            .into_iter()
            .collect();

        assert_eq!(explode_hmm(&hmm, &names, &out).unwrap(), 2);
        assert_eq!(explode_sto(&sto, &names, &out).unwrap(), 2);

        for family in ["alpha", "beta"] {
            assert!(out.join(family).join("query.hmm").is_file(), "{family} hmm");
            assert!(out.join(family).join("query.sto").is_file(), "{family} sto");
        }

        std::fs::remove_dir_all(hmm.parent().unwrap()).ok();
    }

    #[test]
    fn explode_sto_writes_one_directory_per_named_record() {
        let src = tmp("explode.sto", TWO_ALIGNMENTS.as_bytes());
        let out = src.with_file_name("split");

        let names: HashSet<String> = ["beta".to_string()].into_iter().collect();
        assert_eq!(explode_sto(&src, &names, &out).unwrap(), 1);
        assert!(!out.join("alpha").exists());

        let got = ids_and_rows(&out.join("beta").join("query.sto"));
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].0, "beta");
        assert_eq!(got[0].1, vec!["s3 WWWWWWWW"]);

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn subset_sto_keeps_whole_blocks() {
        let src = tmp("subset.sto", TWO_ALIGNMENTS.as_bytes());
        let dst = src.with_file_name("out.sto");

        let names: HashSet<String> = ["alpha".to_string()].into_iter().collect();
        assert_eq!(subset_sto(&src, &names, &dst).unwrap(), 1);

        let got = ids_and_rows(&dst);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].0, "alpha");
        assert_eq!(got[0].1, vec!["s1 ACGTACGT", "s2 ACGTACGT"]);

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }

    #[test]
    fn both_explodes_reject_names_that_escape_the_directory() {
        let hmm = tmp("escape.hmm", two_models().as_bytes());
        let sto = hmm.with_file_name("escape.sto");
        std::fs::write(&sto, TWO_ALIGNMENTS).unwrap();

        let out = hmm.with_file_name("split");
        let names: HashSet<String> = ["../alpha".to_string()].into_iter().collect();

        for err in [
            explode_hmm(&hmm, &names, &out).unwrap_err().to_string(),
            explode_sto(&sto, &names, &out).unwrap_err().to_string(),
        ] {
            assert!(err.contains("cannot be used as a file name"), "got: {err}");
        }

        std::fs::remove_dir_all(hmm.parent().unwrap()).ok();
    }

    /// A latin-1 byte in an author name is what pfam.sto actually holds. The
    /// family still travels, and its alignment is untouched.
    #[test]
    fn an_alignment_that_is_not_utf8_still_travels() {
        let mut body = TWO_ALIGNMENTS.as_bytes().to_vec();
        body.extend_from_slice(
            b"# STOCKHOLM 1.0\n#=GF ID gamma\n#=GF RA Ant\xf4nio RV;\ns4 WWWWAAAA\n//\n",
        );

        let src = tmp("badbyte.sto", &body);
        let dst = src.with_file_name("out.sto");

        let names: HashSet<String> = ["gamma".to_string()].into_iter().collect();
        assert_eq!(subset_sto(&src, &names, &dst).unwrap(), 1);

        let got = ids_and_rows(&dst);
        assert_eq!(got[0].0, "gamma");
        assert_eq!(got[0].1, vec!["s4 WWWWAAAA"], "the alignment is untouched");

        // only the offending byte is replaced
        let text = std::fs::read_to_string(&dst).unwrap();
        assert!(text.contains("Ant?nio RV;"), "got: {text}");

        std::fs::remove_dir_all(src.parent().unwrap()).ok();
    }
}
