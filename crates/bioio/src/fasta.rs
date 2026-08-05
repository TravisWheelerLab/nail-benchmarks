use std::{
    fmt::Display,
    io::{BufRead, BufReader, Read, Seek},
    path::Path,
};

use anyhow::{anyhow, Context};
use indexmap::IndexMap;

#[derive(Default, Clone, PartialEq, Eq, Hash)]
pub struct FastaRecord {
    pub name: String,
    pub extra: String,
    pub seq: String,
}

impl FastaRecord {
    pub fn reverse(&mut self) {
        self.seq = self.seq.chars().rev().collect::<String>();
    }
}

impl Display for FastaRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, ">{} {}", self.name, self.extra)?;
        let mut chunks = self
            .seq
            .as_bytes()
            .chunks(60)
            .map(|c| std::str::from_utf8(c).unwrap());

        if let Some(last) = chunks.next_back() {
            chunks.try_for_each(|c| writeln!(f, "{c}"))?;
            write!(f, "{last}")
        } else {
            Ok(())
        }
    }
}

#[derive(Default)]
pub struct Fasta {
    pub records: IndexMap<String, FastaRecord>,
}

impl Fasta {
    pub fn parse<R: Read>(buf: R) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut fasta = Self::default();
        let mut rec = FastaRecord::default();
        for line in reader.lines() {
            let line = line?;
            if let Some(line) = line.strip_prefix(">") {
                if !rec.name.is_empty() {
                    fasta.records.insert(rec.name.to_string(), rec);
                }
                rec = FastaRecord::default();

                let mut tokens = line.splitn(2, char::is_whitespace);
                rec.name = tokens.next().ok_or(anyhow!("no name"))?.to_string();
                rec.extra = tokens.next().unwrap_or_default().to_string();
            } else {
                rec.seq.push_str(&line)
            }
        }

        if rec != FastaRecord::default() {
            fasta.records.insert(rec.name.to_string(), rec);
        }

        Ok(fasta)
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path.as_ref())?;
        Self::parse(f).with_context(|| format!("failed to parse Fasta from: {:?}", path.as_ref()))
    }
}

pub struct ByteBlock<const L: usize> {
    pub byte_start: u64,
    pub byte_cnts: [u16; L],
}

impl<const L: usize> ByteBlock<L> {
    fn new(byte_start: u64) -> Self {
        Self {
            byte_start,
            byte_cnts: [0; L],
        }
    }
}

pub struct FastaByteIndex<R, const L: usize>
where
    R: Read + Seek,
{
    pub reader: BufReader<R>,
    pub blocks: Vec<ByteBlock<L>>,
    pub size: usize,
}

impl<R, const L: usize> FastaByteIndex<R, L>
where
    R: Read + Seek,
{
    pub fn new(buf: R) -> anyhow::Result<Self> {
        let mut reader = BufReader::new(buf);
        reader.seek(std::io::SeekFrom::Start(0))?;

        let mut total_record_cnt = 0usize;
        let mut block_record_cnt = 0usize;

        let mut record_byte_cnt = 0usize;
        let mut block_byte_start = 0usize;
        let mut byte_offset = 0usize;

        let mut blocks = vec![];
        let mut block: ByteBlock<L> = ByteBlock::new(0);

        for line in (&mut reader).lines() {
            let line = line?;

            if line.starts_with('>') {
                block_byte_start = byte_offset;
                block.byte_cnts[block_record_cnt] = record_byte_cnt as u16;
                block_record_cnt += 1;
                total_record_cnt += 1;
                record_byte_cnt = 0;
            }

            if block_record_cnt == L {
                blocks.push(block);
                block = ByteBlock::new(block_byte_start as u64);
                block_record_cnt = 0;
            }

            record_byte_cnt += line.len() + 1;
            byte_offset += line.len() + 1;
        }

        block.byte_cnts[block_record_cnt] = record_byte_cnt as u16;
        blocks.push(block);

        Ok(Self {
            reader,
            blocks,
            size: total_record_cnt,
        })
    }

    pub fn get_record(&mut self, index: usize) -> anyhow::Result<FastaRecord> {
        let s = self.get(index)?;
        let lines: Vec<&str> = s.lines().collect();

        let header = lines[0];
        let seq = lines
            .iter()
            .skip(1)
            .copied()
            .collect::<Vec<&str>>()
            .join("");

        let (name, extra) = header
            .split_once(char::is_whitespace)
            .unwrap_or((header, ""));

        Ok(FastaRecord {
            name: name[1..].to_string(),
            extra: extra.to_string(),
            seq,
        })
    }

    pub fn get(&mut self, index: usize) -> anyhow::Result<String> {
        let block = &self.blocks[index / L];
        let block_idx = index % L;

        let before = block
            .byte_cnts
            .iter()
            .take(block_idx)
            .map(|&l| l as u64)
            .sum::<u64>();

        let n_read = block.byte_cnts[block_idx] as usize;

        let mut buf = vec![0u8; n_read];
        self.reader
            .seek(std::io::SeekFrom::Start(block.byte_start + before))?;
        self.reader.read_exact(&mut buf)?;

        String::from_utf8(buf).context("failed to produce string")
    }
}

// ---------------------------------------------------------------- shaping

use std::io::{BufWriter, Write};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};

/// Deal a fasta into `n` shards, reshuffling the destination order every `n`
/// records so shards stay comparable in composition rather than reflecting
/// whatever order the source happened to be in.
pub fn split(fa_path: &Path, n: usize, out_dir: &Path, seed: u64) -> anyhow::Result<()> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    let mut writers = Vec::with_capacity(n);
    for i in 1..=n {
        let path = out_dir.join(format!("{i}.fa"));
        let file = std::fs::File::create(&path)
            .with_context(|| format!("failed to create {}", path.display()))?;
        writers.push(BufWriter::new(file));
    }

    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(
        std::fs::File::open(fa_path)
            .with_context(|| format!("failed to open {}", fa_path.display()))?,
    )?;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut order: Vec<usize> = (0..n).collect();

    for i in 1..=index.size {
        let j = i % n;
        if j == 0 {
            order.shuffle(&mut rng);
        }
        let seq = index.get(i)?;
        write!(&mut writers[order[j]], "{seq}")?;
    }

    for mut w in writers {
        w.flush()?;
    }

    Ok(())
}

/// Write a copy of `fa_path` with every sequence reversed. Reversed sequences
/// keep the composition of the original but destroy its homology, which makes
/// them usable as decoys when calibrating score cutoffs.
pub fn reverse(fa_path: &Path, out_path: &Path) -> anyhow::Result<()> {
    if let Some(dir) = out_path.parent() {
        std::fs::create_dir_all(dir)?;
    }

    let mut out = BufWriter::new(
        std::fs::File::create(out_path)
            .with_context(|| format!("failed to create {}", out_path.display()))?,
    );

    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(
        std::fs::File::open(fa_path)
            .with_context(|| format!("failed to open {}", fa_path.display()))?,
    )?;

    for i in 1..=index.size {
        let mut rec = index.get_record(i)?;
        rec.reverse();
        writeln!(out, "{rec}")?;
    }

    out.flush()?;
    Ok(())
}

/// Write the first `n` records of a fasta to `out_path`.
pub fn sample_to(fa_path: &Path, n: usize, out_path: &Path) -> anyhow::Result<()> {
    if let Some(dir) = out_path.parent() {
        std::fs::create_dir_all(dir)?;
    }

    let mut index: FastaByteIndex<_, 64> = FastaByteIndex::new(
        std::fs::File::open(fa_path)
            .with_context(|| format!("failed to open {}", fa_path.display()))?,
    )?;
    let mut out = BufWriter::new(std::fs::File::create(out_path)?);

    for i in 1..=n.min(index.size) {
        write!(out, "{}", index.get(i)?)?;
    }

    out.flush()?;
    Ok(())
}

/// Number of records in a fasta, counted without holding the file in memory.
pub fn count(path: impl AsRef<Path>) -> anyhow::Result<usize> {
    use std::io::BufRead;

    let reader = std::io::BufReader::new(std::fs::File::open(path.as_ref())?);
    let mut n = 0usize;

    for line in reader.lines() {
        if line?.starts_with('>') {
            n += 1;
        }
    }

    Ok(n)
}

/// Residue count of a fasta holding exactly one sequence.
pub fn residue_len(path: impl AsRef<Path>) -> anyhow::Result<usize> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    let mut headers = 0usize;
    let mut len = 0usize;

    for line in text.lines() {
        if line.starts_with('>') {
            headers += 1;
        } else {
            len += line.chars().filter(|c| !c.is_whitespace()).count();
        }
    }

    if headers != 1 {
        anyhow::bail!(
            "expected exactly one sequence in {}, found {headers}",
            path.display()
        );
    }

    Ok(len)
}

/// Generate decoy records by drawing a subsequence from `source` matched to the
/// length of a randomly chosen record in `lengths`, then shuffling it.
///
/// Shuffling preserves amino acid composition while destroying any real
/// homology, which is what makes a decoy a fair negative rather than simply an
/// unrelated sequence.
pub fn decoys(
    source: &Fasta,
    lengths: &[usize],
    count: usize,
    rng: &mut StdRng,
) -> anyhow::Result<Vec<FastaRecord>> {
    if source.records.is_empty() || lengths.is_empty() {
        anyhow::bail!("cannot generate decoys from an empty source");
    }

    let n_src = source.records.len();
    let mut out = Vec::with_capacity(count);
    let mut src_bytes: &[u8] = &[];

    for idx in 0..count {
        let decoy_len = lengths[rng.random_range(0..lengths.len())];

        // keep drawing until a source sequence is long enough to cut from
        while src_bytes.len() < decoy_len {
            src_bytes = source
                .records
                .get_index(rng.random_range(0..n_src))
                .context("bad source index")?
                .1
                .seq
                .as_bytes();
        }

        let start = rng.random_range(0..=src_bytes.len() - decoy_len);
        let mut sample: Vec<u8> = src_bytes[start..start + decoy_len].to_vec();
        sample.shuffle(rng);

        out.push(FastaRecord {
            name: format!("decoy{idx}"),
            extra: String::new(),
            seq: std::str::from_utf8(&sample)
                .context("decoy sequence was not utf8")?
                .to_string(),
        });

        src_bytes = &[];
    }

    Ok(out)
}
