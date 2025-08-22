use std::io::{BufRead, BufReader, Read};

use super::Hit;

// #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
// # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
// #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
// ABATE/15/20-150      -            191 ABATE                PF07336.15   131   2.4e-18   71.7   4.8   1   1   2.1e-22     3e-18   71.3   4.8     3   130    21   149    19   150 0.92 domain: Q7ND30_GLOVI/3-133
// ABATE/21/151-279     -            440 ABATE                PF07336.15   131   2.8e-15   61.8   0.4   1   1     4e-19   5.7e-15   60.7   0.4     3   129   151   279   150   281 0.93 domain: C0Z6E1_BREBN/4-132
// ABATE/20/142-284     -            356 ABATE                PF07336.15   131   2.5e-13   55.4   0.0   1   1   4.9e-17     7e-13   54.0   0.0     1   130   142   283   142   284 0.94 domain: Q5WAM5_ALKCK/16-158
// ABATE/14/178-303     -            376 ABATE                PF07336.15   131   2.8e-13   55.2   7.5   1   1   4.3e-17   6.1e-13   54.2   7.5     4   130   176   302   173   303 0.92 domain: H6NTN6_9BACL/5-130

pub fn parse_domtbl<R: Read>(buf: R) -> anyhow::Result<Vec<Hit>> {
    let reader = BufReader::new(buf);

    let mut hits = vec![];
    for line in reader.lines() {
        let line = line.unwrap_or_default();
        if line.starts_with('#') {
            continue;
        }

        let tokens = line.split_whitespace().collect::<Vec<_>>();

        hits.push(Hit {
            query: tokens[3].to_string(),
            target: tokens[0].to_string(),
            query_start: tokens[15].parse()?,
            query_end: tokens[16].parse()?,
            target_start: tokens[17].parse()?,
            target_end: tokens[18].parse()?,
            score: tokens[13].parse()?,
            e_value: tokens[12].parse()?,
        })
    }

    Ok(hits)
}
