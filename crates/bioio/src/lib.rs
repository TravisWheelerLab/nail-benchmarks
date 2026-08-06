pub mod aggregate;
pub mod fasta;
pub mod hmm;
pub mod mmseqs;
pub mod split;
pub mod stockholm;
pub mod tbl;

/// Reject record identifiers that cannot safely become file names.
///
/// The `explode` functions turn identifiers into paths, and identifiers come
/// out of whatever file was handed in. A separator or `..` would put the
/// output somewhere the caller did not ask for.
pub(crate) fn check_file_names(names: &std::collections::HashSet<String>) -> anyhow::Result<()> {
    for name in names {
        if name.is_empty()
            || name == ".."
            || name.contains(std::path::MAIN_SEPARATOR)
            || name.contains('/')
        {
            anyhow::bail!("record identifier {name:?} cannot be used as a file name");
        }
    }
    Ok(())
}
