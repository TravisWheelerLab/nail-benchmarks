use std::io::{BufRead, BufReader, Read};

pub struct Hit {
    pub query: String,
    pub target: String,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub score: f64,
    pub e_value: f64,
}

// #                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
// # target name  accession   tlen query name  accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
// #------------- ---------- ----- ----------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
pub fn parse_hmmer_domtbl<R: Read>(buf: R) -> anyhow::Result<Vec<Hit>> {
    parse_tbl::<R, 3, 0, 15, 16, 17, 18, 13, 12>(buf)
}

// #                target target query query       comp        cell
// # target  query  start  end    start end   score bias evalue frac
// # ------- ------ ------ ------ ----- ----- ----- ---- ------ -----
pub fn parse_nail_tbl<R: Read>(buf: R) -> anyhow::Result<Vec<Hit>> {
    parse_tbl::<R, 1, 0, 4, 5, 2, 3, 6, 8>(buf)
}

// # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
pub fn parse_blast_tbl<R: Read>(buf: R) -> anyhow::Result<Vec<Hit>> {
    parse_tbl::<R, 0, 1, 6, 7, 8, 9, 11, 10>(buf)
}

pub fn parse_tbl<
    R: Read,
    const QUERY: usize,
    const TARGET: usize,
    const Q_START: usize,
    const Q_END: usize,
    const T_START: usize,
    const T_END: usize,
    const SCORE: usize,
    const E_VALUE: usize,
>(
    buf: R,
) -> anyhow::Result<Vec<Hit>> {
    let reader = BufReader::new(buf);

    let mut hits = vec![];
    for line in reader.lines() {
        let line = line.unwrap_or_default();
        if line.starts_with('#') {
            continue;
        }

        let tokens = line.split_whitespace().collect::<Vec<_>>();

        hits.push(Hit {
            query: tokens[QUERY].to_string(),
            target: tokens[TARGET].to_string(),
            query_start: tokens[Q_START].parse()?,
            query_end: tokens[Q_END].parse()?,
            target_start: tokens[T_START].parse()?,
            target_end: tokens[T_END].parse()?,
            score: tokens[SCORE].parse()?,
            e_value: tokens[E_VALUE].parse()?,
        })
    }

    Ok(hits)
}
