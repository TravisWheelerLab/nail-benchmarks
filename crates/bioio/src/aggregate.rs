//! Random access and sampling across a directory of fasta files treated as one
//! collection.
//!
//! The collection this was built for is ~500GB across 25 files, so nothing here
//! reads more than it has to. Construction only lists files; indexing is
//! explicit and can be persisted, and sampling reads only the records it draws.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

use anyhow::{bail, Context, Result};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};

/// Records per index block. A block stores one byte offset, so this trades
/// index size against how far a lookup scans: at 64, an index costs 0.125
/// bytes per record and a lookup scans at most 64 records.
pub const BLOCK: usize = 64;

const MAGIC: &[u8; 7] = b"BIOIDX\0";
const VERSION: u16 = 1;

/// Bytes read at a time when scanning. Records are small relative to this, so
/// a lookup usually needs one read.
const BUF: usize = 1 << 20;

/// What a file looked like when it was indexed. Any change invalidates the
/// index rather than silently sampling from stale offsets.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Stamp {
    name: String,
    size: u64,
    mtime_ns: u128,
}

impl Stamp {
    fn of(path: &Path) -> Result<Self> {
        let meta = std::fs::metadata(path)
            .with_context(|| format!("failed to stat {}", path.display()))?;
        let mtime_ns = meta
            .modified()?
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);

        Ok(Stamp {
            name: path
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default(),
            size: meta.len(),
            mtime_ns,
        })
    }
}

/// A sparse index over one file: the byte offset of the first record of every
/// block of `BLOCK` records. Blocks always begin on a record boundary.
#[derive(Clone, Debug)]
struct FileIndex {
    stamp: Stamp,
    block_starts: Vec<u64>,
    records: u64,
}

/// A set of fasta files addressed as a single collection.
#[derive(Debug)]
pub struct AggregateFasta {
    paths: Vec<PathBuf>,
    index: Index,
}

/// The built index: per-file sparse offsets plus the cumulative record counts
/// that map a global record number onto a file.
#[derive(Clone, Debug)]
struct Index {
    files: Vec<FileIndex>,
    /// `starts[i]` is the global number of the first record in file `i`, and
    /// `starts[n]` is the total. This is the virtual index: ranges, not a
    /// per-record map.
    starts: Vec<u64>,
}

impl Index {
    fn total(&self) -> u64 {
        *self.starts.last().unwrap_or(&0)
    }

    /// Which file a global record number falls in, and its offset within it.
    fn locate(&self, global: u64) -> Option<(usize, u64)> {
        if global >= self.total() {
            return None;
        }
        // starts is sorted, so the owning file is the last one starting at or
        // before `global`
        let file = match self.starts.binary_search(&global) {
            Ok(i) => i,
            Err(i) => i - 1,
        };
        Some((file, global - self.starts[file]))
    }
}

/// Extensions treated as fasta when none are given.
pub const DEFAULT_EXTENSIONS: [&str; 2] = ["fa", "fasta"];

/// Assembles an [`AggregateFasta`].
///
/// Sources accumulate: `dir` and `path` may each be called any number of times
/// and mixed, so a collection can span directories. Duplicates are removed, so
/// naming a file both directly and via its directory is harmless.
pub struct AggregateFastaBuilder {
    dirs: Vec<PathBuf>,
    paths: Vec<PathBuf>,
    extensions: Vec<String>,
    index_path: Option<PathBuf>,
    allow_overwrite: bool,
}

impl AggregateFastaBuilder {
    fn new() -> Self {
        AggregateFastaBuilder {
            dirs: Vec::new(),
            paths: Vec::new(),
            extensions: DEFAULT_EXTENSIONS.iter().map(|s| s.to_string()).collect(),
            index_path: None,
            allow_overwrite: false,
        }
    }

    /// Add every fasta in a directory. Not recursive.
    pub fn dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.dirs.push(dir.into());
        self
    }

    /// Add one file, whatever its extension.
    pub fn path(mut self, path: impl Into<PathBuf>) -> Self {
        self.paths.push(path.into());
        self
    }

    /// Add extensions that count as fasta when scanning a directory.
    ///
    /// Additive: `fa` and `fasta` are always recognised. A leading dot is
    /// optional, and everything is lowercased before comparison, so `pep`,
    /// `.pep` and `.PEP` are one extension rather than three.
    pub fn extensions<I, S>(mut self, exts: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.extensions.extend(
            exts.into_iter()
                .map(|e| e.as_ref().trim_start_matches('.').to_ascii_lowercase()),
        );
        self.extensions.sort();
        self.extensions.dedup();
        self
    }

    /// Read the index from `path` if it is there and still matches the sources,
    /// otherwise build it and write it there.
    ///
    /// Omit this and the index is built in memory and never written.
    pub fn index(mut self, path: impl Into<PathBuf>) -> Self {
        self.index_path = Some(path.into());
        self
    }

    /// Rebuild and replace an index that no longer matches its sources.
    /// Without this, a mismatch is an error.
    pub fn allow_overwrite(mut self) -> Self {
        self.allow_overwrite = true;
        self
    }

    /// Extensions currently recognised, sorted and deduplicated.
    pub fn extensions_list(&self) -> &[String] {
        &self.extensions
    }

    /// Resolve the sources and produce an indexed collection.
    pub fn build(self) -> Result<AggregateFasta> {
        let paths = self.collect_paths()?;

        let index = match &self.index_path {
            None => build(index_all(&paths)?),
            Some(cache) if cache.exists() => {
                let stored = read_index(cache)
                    .with_context(|| format!("failed to read index {}", cache.display()))?;

                if stamps_match(&stored, &paths)? {
                    stored
                } else if self.allow_overwrite {
                    eprintln!(
                        "warning: index {} no longer matches its sources; rebuilding it",
                        cache.display()
                    );
                    let fresh = build(index_all(&paths)?);
                    write_index(cache, &fresh)?;
                    fresh
                } else {
                    bail!(
                        "index {} no longer matches its sources; \
                         pass allow_overwrite to rebuild it",
                        cache.display()
                    );
                }
            }
            Some(cache) => {
                let fresh = build(index_all(&paths)?);
                write_index(cache, &fresh)?;
                fresh
            }
        };

        Ok(AggregateFasta { paths, index })
    }

    /// Expand directories, keep explicit paths, then sort and deduplicate.
    fn collect_paths(&self) -> Result<Vec<PathBuf>> {
        let mut out: Vec<PathBuf> = Vec::new();

        for dir in &self.dirs {
            let entries = std::fs::read_dir(dir)
                .with_context(|| format!("failed to read directory {}", dir.display()))?;

            for entry in entries {
                let path = entry?.path();
                let matches = path.extension().is_some_and(|e| {
                    let e = e.to_string_lossy().to_ascii_lowercase();
                    self.extensions.iter().any(|want| *want == e)
                });

                if path.is_file() && matches {
                    out.push(path);
                }
            }
        }

        for path in &self.paths {
            if !path.is_file() {
                bail!("{} is not a readable file", path.display());
            }
            out.push(path.clone());
        }

        // canonicalise before deduplicating, so the same file reached by two
        // routes is not counted twice, which would double its records in the
        // global numbering
        let mut canonical: Vec<PathBuf> = out
            .iter()
            .map(|p| {
                p.canonicalize()
                    .with_context(|| format!("failed to resolve {}", p.display()))
            })
            .collect::<Result<_>>()?;

        // global record numbering depends on a stable order
        canonical.sort();
        canonical.dedup();

        if canonical.is_empty() {
            bail!("no fasta files found; add a dir() or a path()");
        }

        Ok(canonical)
    }
}

impl AggregateFasta {
    pub fn builder() -> AggregateFastaBuilder {
        AggregateFastaBuilder::new()
    }

    pub fn files(&self) -> &[PathBuf] {
        &self.paths
    }

    /// Total records across the collection.
    pub fn len(&self) -> u64 {
        self.index.total()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Draw `n` records uniformly at random without replacement and write them
    /// verbatim to `out`.
    ///
    /// Returns how many were written, which is fewer than `n` only when the
    /// collection is smaller than the request.
    pub fn sample(&self, n: usize, seed: u64, out: &mut impl Write) -> Result<usize> {
        let index = &self.index;
        let total = index.total();
        let n = if n as u64 > total {
            eprintln!(
                "warning: asked for {n} records but the collection holds {total}; using all of them"
            );
            total as usize
        } else {
            n
        };

        let mut rng = StdRng::seed_from_u64(seed);

        // Draw without replacement, in random order.
        //
        // The order is deliberately not sorted: the emitted records must be a
        // random permutation, not the collection's own order with gaps. Reads
        // are scattered either way — draws over a large collection land far
        // apart no matter how they are ordered — so there is little to reclaim
        // by sorting, and doing so would make the output non-random.
        //
        // Two strategies, because materialising 0..total is untenable at scale
        // — a 2e9-record collection would want a 20GB vector to draw a few
        // thousand records.
        let draw: Vec<u64> = if (n as u64).saturating_mul(3) >= total {
            // taking a large fraction: shuffling everything is cheapest
            let mut all: Vec<u64> = (0..total).collect();
            all.shuffle(&mut rng);
            all.truncate(n);
            all
        } else {
            // taking a sparse subset: reject duplicates, ~1.5n draws expected
            let mut seen = std::collections::HashSet::with_capacity(n);
            while seen.len() < n {
                seen.insert(rng.random_range(0..total));
            }
            // hash order is not reproducible, so impose the seeded one
            let mut v: Vec<u64> = seen.into_iter().collect();
            v.sort_unstable();
            v.shuffle(&mut rng);
            v
        };

        // one reader per file, opened on first use and kept: draws arrive in
        // random order, so reopening on every file change would thrash
        let mut readers: Vec<Option<BufReader<File>>> = (0..self.paths.len()).map(|_| None).collect();
        let mut written = 0usize;

        for global in draw {
            let (file_idx, local) = index
                .locate(global)
                .context("drew a record outside the collection")?;

            if readers[file_idx].is_none() {
                let f = File::open(&self.paths[file_idx])
                    .with_context(|| format!("failed to open {}", self.paths[file_idx].display()))?;
                readers[file_idx] = Some(BufReader::new(f));
            }

            let reader = readers[file_idx].as_mut().expect("reader was just opened");
            let bytes = read_record(reader, &index.files[file_idx], local)?;
            out.write_all(&bytes)?;
            written += 1;
        }

        Ok(written)
    }
}

/// Index one file, recording the byte offset of every `BLOCK`th record.
fn index_file(path: &Path) -> Result<FileIndex> {
    let stamp = Stamp::of(path)?;
    let mut reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );

    let mut buf = vec![0u8; BUF];
    let mut block_starts = Vec::new();
    let mut records: u64 = 0;
    let mut pos: u64 = 0;
    // whether the previous buffer ended on a newline, so a '>' opening the next
    // one still counts as a record start
    let mut at_line_start = true;

    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }

        // a record opens at a '>' that begins a line, so jump between newlines
        // rather than inspecting every byte
        if at_line_start && buf[0] == b'>' {
            if records % BLOCK as u64 == 0 {
                block_starts.push(pos);
            }
            records += 1;
        }

        for i in memchr::memchr_iter(b'\n', &buf[..n]) {
            if i + 1 < n && buf[i + 1] == b'>' {
                if records % BLOCK as u64 == 0 {
                    block_starts.push(pos + i as u64 + 1);
                }
                records += 1;
            }
        }

        at_line_start = buf[n - 1] == b'\n';
        pos += n as u64;
    }

    Ok(FileIndex {
        stamp,
        block_starts,
        records,
    })
}

/// Index every file, concurrently.
///
/// Uses rayon's global pool, so a binary that sizes it from a `--threads` flag
/// controls this too. Indexing is disk-bound: on a 528GB collection throughput
/// went 3.1 -> 5.8 -> 8.5 GB/s at 1, 2 and 4 readers, then stayed flat through
/// 20, so there is nothing to gain past a handful.
fn index_all(paths: &[PathBuf]) -> Result<Vec<FileIndex>> {
    use rayon::prelude::*;

    // collect preserves order, which global record numbering depends on
    paths.par_iter().map(|p| index_file(p)).collect()
}

/// Whether a stored index still describes these files unchanged.
fn stamps_match(index: &Index, paths: &[PathBuf]) -> Result<bool> {
    if index.files.len() != paths.len() {
        return Ok(false);
    }

    for (stored, path) in index.files.iter().zip(paths) {
        if stored.stamp != Stamp::of(path)? {
            return Ok(false);
        }
    }

    Ok(true)
}

fn build(files: Vec<FileIndex>) -> Index {
    let mut starts = Vec::with_capacity(files.len() + 1);
    let mut acc = 0u64;
    for f in &files {
        starts.push(acc);
        acc += f.records;
    }
    starts.push(acc);

    Index { files, starts }
}

/// Read record `local` out of an indexed file, verbatim.
///
/// Blocks begin and end on record boundaries, so the next block's start is
/// exactly where this block ends. Reading that span gets the whole block in one
/// go — about 17KB at the default block length — rather than scanning forward
/// through arbitrarily sized buffers.
fn read_record(reader: &mut BufReader<File>, index: &FileIndex, local: u64) -> Result<Vec<u8>> {
    let block = (local / BLOCK as u64) as usize;
    let within = (local % BLOCK as u64) as usize;

    let start = *index
        .block_starts
        .get(block)
        .context("record outside the indexed range")?;

    reader.seek(SeekFrom::Start(start))?;

    // read exactly this block, or to the end of the file for the last one
    let mut buf = Vec::new();
    match index.block_starts.get(block + 1) {
        Some(&next) => {
            buf.resize((next - start) as usize, 0);
            reader.read_exact(&mut buf)?;
        }
        None => {
            reader.read_to_end(&mut buf)?;
        }
    }

    // record starts within the block: the first byte, then every '>' opening a
    // line after it
    let mut starts = Vec::with_capacity(BLOCK);
    starts.push(0usize);
    for i in memchr::memchr_iter(b'\n', &buf) {
        if i + 1 < buf.len() && buf[i + 1] == b'>' {
            starts.push(i + 1);
        }
    }

    let begin = *starts
        .get(within)
        .context("record not found in its block; the index may not match the file")?;
    let end = starts.get(within + 1).copied().unwrap_or(buf.len());

    buf.truncate(end);
    buf.drain(..begin);
    Ok(buf)
}

// ------------------------------------------------------------- persistence

fn write_index(path: &Path, index: &Index) -> Result<()> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let mut out = std::io::BufWriter::new(File::create(path)?);

    out.write_all(MAGIC)?;
    out.write_all(&VERSION.to_le_bytes())?;
    out.write_all(&(BLOCK as u32).to_le_bytes())?;
    out.write_all(&(index.files.len() as u32).to_le_bytes())?;

    for f in &index.files {
        let name = f.stamp.name.as_bytes();
        out.write_all(&(name.len() as u16).to_le_bytes())?;
        out.write_all(name)?;
        out.write_all(&f.stamp.size.to_le_bytes())?;
        out.write_all(&f.stamp.mtime_ns.to_le_bytes())?;
        out.write_all(&f.records.to_le_bytes())?;
        out.write_all(&(f.block_starts.len() as u64).to_le_bytes())?;
        for b in &f.block_starts {
            out.write_all(&b.to_le_bytes())?;
        }
    }

    Ok(())
}

fn read_index(path: &Path) -> Result<Index> {
    let mut r = BufReader::new(File::open(path)?);

    let mut magic = [0u8; 7];
    r.read_exact(&mut magic)?;
    if &magic != MAGIC {
        bail!("{} is not an index", path.display());
    }

    if read_u16(&mut r)? != VERSION {
        bail!("{} was written by a different version", path.display());
    }
    if read_u32(&mut r)? as usize != BLOCK {
        bail!("{} uses a different block length", path.display());
    }

    let n_files = read_u32(&mut r)? as usize;
    let mut files = Vec::with_capacity(n_files);

    for _ in 0..n_files {
        let name_len = read_u16(&mut r)? as usize;
        let mut name = vec![0u8; name_len];
        r.read_exact(&mut name)?;

        let stamp = Stamp {
            name: String::from_utf8(name).context("index holds a non-utf8 file name")?,
            size: read_u64(&mut r)?,
            mtime_ns: read_u128(&mut r)?,
        };

        let records = read_u64(&mut r)?;
        let n_blocks = read_u64(&mut r)? as usize;
        let mut block_starts = Vec::with_capacity(n_blocks);
        for _ in 0..n_blocks {
            block_starts.push(read_u64(&mut r)?);
        }

        files.push(FileIndex {
            stamp,
            block_starts,
            records,
        });
    }

    Ok(build(files))
}

fn read_u16(r: &mut impl Read) -> Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}

fn read_u32(r: &mut impl Read) -> Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64(r: &mut impl Read) -> Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn read_u128(r: &mut impl Read) -> Result<u128> {
    let mut b = [0u8; 16];
    r.read_exact(&mut b)?;
    Ok(u128::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A collection whose records are individually identifiable, so a sample
    /// can be checked for correctness rather than just for size.
    fn collection(name: &str, files: usize, per_file: usize) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("bioio-agg-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        for f in 0..files {
            let mut body = String::new();
            for r in 0..per_file {
                // varying sequence length, so nothing can rely on fixed strides
                let len = 5 + (f * 7 + r * 13) % 40;
                body.push_str(&format!(">f{f}r{r}\n{}\n", "A".repeat(len)));
            }
            std::fs::write(dir.join(format!("{f:02}.fa")), body).unwrap();
        }

        dir
    }

    fn names(text: &str) -> Vec<String> {
        text.lines()
            .filter(|l| l.starts_with('>'))
            .map(|l| l[1..].to_string())
            .collect()
    }

    #[test]
    fn indexes_every_record_across_files() {
        let dir = collection("count", 3, 150);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        assert_eq!(agg.files().len(), 3);
        assert_eq!(agg.len(), 450);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn sample_draws_whole_records_without_replacement() {
        let dir = collection("draw", 3, 100);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        let mut out = Vec::new();
        let n = agg.sample(120, 42, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let got = names(&text);

        assert_eq!(n, 120);
        assert_eq!(got.len(), 120, "every record needs its header");

        let unique: std::collections::HashSet<_> = got.iter().collect();
        assert_eq!(unique.len(), 120, "sampling must be without replacement");

        // sequences must arrive intact, not truncated at a block boundary
        for line in text.lines().filter(|l| !l.starts_with('>')) {
            assert!(line.len() >= 5 && line.chars().all(|c| c == 'A'), "bad seq: {line:?}");
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn sample_reaches_every_file_and_both_ends() {
        let dir = collection("spread", 4, 100);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        let mut out = Vec::new();
        agg.sample(400, 7, &mut out).unwrap();
        let got = names(&String::from_utf8(out).unwrap());

        // drawing the whole collection must yield exactly it
        assert_eq!(got.len(), 400);
        for f in 0..4 {
            assert!(got.contains(&format!("f{f}r0")), "missing first record of file {f}");
            assert!(got.contains(&format!("f{f}r99")), "missing last record of file {f}");
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn records_spanning_block_boundaries_read_intact() {
        // more records than one block, so lookups must scan within a block
        let dir = collection("blocks", 1, BLOCK * 3 + 7);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        let mut out = Vec::new();
        agg.sample(BLOCK * 3 + 7, 1, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        // the file's own content, reordered, is exactly what should come back
        let source = std::fs::read_to_string(dir.join("00.fa")).unwrap();
        let mut want = names(&source);
        let mut got = names(&text);
        want.sort();
        got.sort();
        assert_eq!(got, want);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn asking_for_too_many_clamps() {
        let dir = collection("clamp", 2, 10);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        let mut out = Vec::new();
        let n = agg.sample(500, 3, &mut out).unwrap();
        assert_eq!(n, 20);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_index_file_round_trips() {
        let dir = collection("cache", 2, 80);
        let idx = dir.join(".index");

        let agg = AggregateFasta::builder().dir(&dir).index(&idx).build().unwrap();
        assert!(idx.exists(), "index() should have written the file");
        assert_eq!(agg.len(), 160);

        // reloading must reuse the file rather than rebuild
        let again = AggregateFasta::builder().dir(&dir).index(&idx).build().unwrap();
        assert_eq!(again.len(), 160);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_stale_index_errors_unless_overwrite_is_allowed() {
        let dir = collection("stale", 2, 80);
        let idx = dir.join(".index");

        AggregateFasta::builder().dir(&dir).index(&idx).build().unwrap();

        // changing a source leaves the stored stamp behind
        std::thread::sleep(std::time::Duration::from_millis(10));
        std::fs::write(dir.join("00.fa"), ">new\nAAAA\n").unwrap();

        let err = AggregateFasta::builder()
            .dir(&dir)
            .index(&idx)
            .build()
            .unwrap_err()
            .to_string();
        assert!(err.contains("no longer matches"), "unexpected: {err}");

        let rebuilt = AggregateFasta::builder()
            .dir(&dir)
            .index(&idx)
            .allow_overwrite()
            .build()
            .unwrap();
        assert_eq!(rebuilt.len(), 81, "should have rebuilt from the changed sources");

        // and the replacement must load cleanly without the flag
        let reloaded = AggregateFasta::builder().dir(&dir).index(&idx).build().unwrap();
        assert_eq!(reloaded.len(), 81);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn sources_accumulate_and_deduplicate() {
        let dir = collection("sources", 3, 10);
        let extra = dir.join("other");
        std::fs::create_dir_all(&extra).unwrap();
        std::fs::write(extra.join("x.fasta"), ">x0\nAAAA\n>x1\nCCCC\n").unwrap();

        // the same file by two routes, plus a second directory
        let agg = AggregateFasta::builder()
            .dir(&dir)
            .path(dir.join("00.fa"))
            .dir(&extra)
            .build()
            .unwrap();

        assert_eq!(agg.files().len(), 4, "00.fa must not be counted twice");
        assert_eq!(agg.len(), 32, "3 x 10 records plus the 2 in x.fasta");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn extensions_append_to_the_defaults_and_ignore_case() {
        let dir = collection("exts", 1, 5);
        std::fs::write(dir.join("upper.FA"), ">u0\nAAAA\n").unwrap();
        std::fs::write(dir.join("odd.pep"), ">p0\nAAAA\n>p1\nCCCC\n").unwrap();

        // .FA counts by default, .pep does not
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();
        assert_eq!(agg.len(), 6);

        // adding .pep keeps the defaults rather than replacing them
        let agg = AggregateFasta::builder()
            .dir(&dir)
            .extensions([".pep"])
            .build()
            .unwrap();
        assert_eq!(agg.len(), 8, "the 6 default-matched plus the 2 in .pep");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn extension_spellings_collapse_to_one() {
        let agg = AggregateFasta::builder()
            .extensions([".PEP", "pep"])
            .extensions(["Pep", "fa", ".FASTA"]);

        // dots, case and repeats across calls all fold together
        assert_eq!(agg.extensions_list(), ["fa", "fasta", "pep"]);
    }

    #[test]
    fn a_missing_path_is_an_error() {
        let dir = collection("missing", 1, 4);

        let err = AggregateFasta::builder()
            .path(dir.join("nope.fa"))
            .build()
            .unwrap_err()
            .to_string();
        assert!(err.contains("not a readable file"), "unexpected: {err}");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn building_with_no_sources_is_an_error() {
        let err = AggregateFasta::builder().build().unwrap_err().to_string();
        assert!(err.contains("no fasta files"), "unexpected: {err}");
    }

    #[test]
    fn output_order_is_a_random_permutation_not_collection_order() {
        let dir = collection("order", 4, 200);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        let mut out = Vec::new();
        agg.sample(200, 5, &mut out).unwrap();
        let got = names(&String::from_utf8(out).unwrap());

        // sorted draws would emit every record of file 0 before any of file 1,
        // and each file's records in increasing order
        let file_of = |n: &String| n[1..n.find('r').unwrap()].parse::<usize>().unwrap();
        let files: Vec<usize> = got.iter().map(file_of).collect();
        assert!(
            files.windows(2).any(|w| w[0] > w[1]),
            "files appear in collection order, so the output was not shuffled"
        );

        // and within one file, record numbers must not be ascending either
        let rec_of = |n: &String| n[n.find('r').unwrap() + 1..].parse::<usize>().unwrap();
        let first: Vec<usize> = got.iter().filter(|n| file_of(n) == 0).map(rec_of).collect();
        assert!(
            first.windows(2).any(|w| w[0] > w[1]),
            "records within a file are ascending, so the output was not shuffled"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn sampling_is_reproducible_for_a_seed() {
        let dir = collection("seed", 2, 50);
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

        let mut a = Vec::new();
        let mut b = Vec::new();
        agg.sample(30, 99, &mut a).unwrap();
        agg.sample(30, 99, &mut b).unwrap();
        assert_eq!(a, b);

        let mut c = Vec::new();
        agg.sample(30, 100, &mut c).unwrap();
        assert_ne!(a, c);

        std::fs::remove_dir_all(&dir).ok();
    }
}
