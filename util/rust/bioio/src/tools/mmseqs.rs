use std::io::{BufRead, BufReader, Read};

use super::Hit;

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
            e_value: tokens[7].parse()?,
        })
    }

    Ok(hits)
}
