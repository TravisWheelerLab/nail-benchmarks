use std::io::{BufRead, BufReader, Read};

use super::Hit;

// #                                                               target target query query       comp        cell
// # target                              query                     start  end    start end   score bias evalue frac
// # ----------------------------------- ------------------------- ------ ------ ----- ----- ----- ---- ------ -----
// DUF4880:327-369|49%|Q6NA38_RHOPA/9-51 DUF4880|Q6NA38_RHOPA/9-51 328    368    2     42    30.2  0.7  2.5e-5 0.043
// DUF4880:37-78|45%|A9BW91_DELAS/19-59  DUF4880|Q6NA38_RHOPA/9-51 38     75     2     39    19.4  1.2  8.6e-2 0.165
// DUF4880:12-54|44%|Q883J9_PSESM/16-58  DUF4880|Q6NA38_RHOPA/9-51 13     50     2     39    17.7  0.0  3.1e-1 0.333
// decoy1239421                          DUF4880|Q6NA38_RHOPA/9-51 83     101    3     21    13.5  0.0  7.3e0  0.155
// decoy1066413                          DUF4880|Q6NA38_RHOPA/9-51 275    295    8     28    12.3  0.0  1.8e1  0.058
// decoy1071907                          DUF4880|Q6NA38_RHOPA/9-51 160    184    10    34    12.1  0.0  2.2e1  0.083
// decoy397873                           DUF4880|Q6NA38_RHOPA/9-51 33     53     16    36    11.7  0.1  2.8e1  0.059

pub fn parse_tbl<R: Read>(buf: R) -> anyhow::Result<Vec<Hit>> {
    let reader = BufReader::new(buf);

    let mut hits = vec![];
    for line in reader.lines() {
        let line = line.unwrap_or_default();
        if line.starts_with('#') {
            continue;
        }

        let tokens = line.split_whitespace().collect::<Vec<_>>();

        hits.push(Hit {
            query: tokens[1].to_string(),
            target: tokens[0].to_string(),
            query_start: tokens[4].parse()?,
            query_end: tokens[5].parse()?,
            target_start: tokens[2].parse()?,
            target_end: tokens[3].parse()?,
            score: tokens[6].parse()?,
            e_value: tokens[8].parse()?,
        })
    }

    Ok(hits)
}
