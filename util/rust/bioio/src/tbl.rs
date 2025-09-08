use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use anyhow::Context;

pub struct Hit {
    pub query: String,
    pub target: String,
    pub query_start: Option<usize>,
    pub query_end: Option<usize>,
    pub target_start: Option<usize>,
    pub target_end: Option<usize>,
    pub score: f64,
    pub e_value: f64,
}

pub trait HitColumns {
    const QUERY: usize;
    const TARGET: usize;
    const Q_START: Option<usize>;
    const Q_END: Option<usize>;
    const T_START: Option<usize>;
    const T_END: Option<usize>;
    const SCORE: usize;
    const E_VALUE: usize;
}

// |      0      |    1     |     2      |     3    |    4   |  5   |  6  |    7    |   8  |  9  |  10 | 11| 12| 13| 14| 15| 16| 17|         18          |
// #                                                                             --- full sequence ---- --- best 1 domain ---- --- domain number estimation
// # target name   accession   query name  accession  E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
// #------------- ---------- ------------ ---------- -------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
pub struct HmmerTable {}
impl HitColumns for HmmerTable {
    const QUERY: usize = 2;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = None;
    const Q_END: Option<usize> = None;
    const T_START: Option<usize> = None;
    const T_END: Option<usize> = None;
    const SCORE: usize = 5;
    const E_VALUE: usize = 4;
}

// |      0      |    1     |  2  |     3     |    4     |  5  |    6    |   7  |  8  | 9 | 10|   11    |    12   |  13  |  14 |  15 |  16 |  17 |  18 |  19 | 20  | 21 |         22          |
// #                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
// # target name  accession   tlen query name  accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
// #------------- ---------- ----- ----------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
pub struct HmmerDomTable {}
impl HitColumns for HmmerDomTable {
    const QUERY: usize = 3;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = Some(15);
    const Q_END: Option<usize> = Some(16);
    const T_START: Option<usize> = Some(17);
    const T_END: Option<usize> = Some(18);
    const SCORE: usize = 13;
    const E_VALUE: usize = 12;
}

//  |   0   |   1  |  2   |   3  |  4  |  5  |  6  |  7 |   8  |  9  |
// #                target target query query       comp        cell
// # target  query  start  end    start end   score bias evalue frac
// # ------- ------ ------ ------ ----- ----- ----- ---- ------ -----
pub struct NailTable {}
impl HitColumns for NailTable {
    const QUERY: usize = 1;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = Some(4);
    const Q_END: Option<usize> = Some(5);
    const T_START: Option<usize> = Some(2);
    const T_END: Option<usize> = Some(3);
    const SCORE: usize = 6;
    const E_VALUE: usize = 8;
}

//          |      0      |        1       |     2     |       3         |     4     |     5    |    6    |   7   |    8    |   9   |   10  |    11    |
// # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
pub struct BlastTable {}
impl HitColumns for BlastTable {
    const QUERY: usize = 0;
    const TARGET: usize = 1;
    const Q_START: Option<usize> = Some(6);
    const Q_END: Option<usize> = Some(7);
    const T_START: Option<usize> = Some(8);
    const T_END: Option<usize> = Some(9);
    const SCORE: usize = 11;
    const E_VALUE: usize = 10;
}

pub struct HitTable {
    pub name: String,
    pub hits: Vec<Hit>,
}

impl HitTable {
    pub fn from_path<P: AsRef<Path>, C: HitColumns>(path: P) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("invalid path")?;

        Self::parse::<File, C>(file, name)
    }

    pub fn parse<R: Read, C: HitColumns>(buf: R, name: &str) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut hits = vec![];
        for line in reader.lines() {
            let line = line.unwrap_or_default();
            if line.starts_with('#') {
                continue;
            }

            let tokens = line.split_whitespace().collect::<Vec<_>>();

            hits.push(Hit {
                query: tokens[C::QUERY].to_string(),
                target: tokens[C::TARGET].to_string(),
                query_start: C::Q_START.map(|i| tokens[i].parse::<usize>()).transpose()?,
                query_end: C::Q_END.map(|i| tokens[i].parse::<usize>()).transpose()?,
                target_start: C::T_START.map(|i| tokens[i].parse::<usize>()).transpose()?,
                target_end: C::T_END.map(|i| tokens[i].parse::<usize>()).transpose()?,
                score: tokens[C::SCORE].parse()?,
                e_value: tokens[C::E_VALUE].parse()?,
            })
        }

        Ok(Self {
            name: name.to_string(),
            hits,
        })
    }
}
