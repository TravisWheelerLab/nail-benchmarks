//! What the benchmarks in this repo agree on.
//!
//! Each benchmark asks its own question and keeps its own analyses. Two things
//! are the same for all of them. The first is the shape of the record a run
//! leaves behind: a [`pail::Table`] sink writes `manifest.tbl`, `parse` reads it
//! back through [`manifest`], and whatever it works out gets written through
//! [`tbl`]. The second is where the programs being benchmarked and the sequence
//! data they run on were put, which is [`tools`].
//!
//! Two smaller things are shared for want of a second home: [`split`], which
//! cuts a query set up for a batch of jobs, and [`nail`], which reads the one
//! column of nail's table that libsail's layout does not carry. [`time`] is
//! for the other way a run can be recorded: timed by a shell rather than by a
//! pipeline, on a machine this workspace never sees.
//!
//! Nothing here decides what counts as a *true* hit. That is the part every
//! benchmark answers differently, and it stays with the benchmark.

pub mod manifest;
pub mod nail;
pub mod split;
pub mod tbl;
pub mod time;
pub mod tools;

/// Replace every byte that is not part of a valid UTF-8 sequence with `?`.
///
/// Pfam's SEED alignments carry author names in Latin-1 -- `pfam.sto` holds an
/// `ô` as a single byte in a `#=GF RA` line -- and libsail wants a record to be
/// text. The bytes are always in `#=GF` metadata rather than in an alignment
/// row or a sequence, so replacing them costs nothing any benchmark reads.
pub fn repair_utf8(bytes: &mut [u8]) {
    let mut from = 0;

    while let Err(e) = std::str::from_utf8(&bytes[from..]) {
        let bad = from + e.valid_up_to();

        // None is a sequence cut off by the end of the input, so there is
        // nothing past it to resume from
        let len = e.error_len().unwrap_or(bytes.len() - bad);

        bytes[bad..bad + len].fill(b'?');
        from = bad + len;
    }
}

/// A profile of `leng` nodes over a two-symbol alphabet, which is the smallest
/// thing libsail's parser accepts.
///
/// Public because the tests that need a profile libsail will read are in two
/// crates: [`split`]'s here, and mgy's around the Pfam cuts.
pub fn profile(name: &str, leng: usize) -> String {
    let mut out = format!(
        "HMMER3/f\nNAME  {name}\nLENG  {leng}\nALPH  amino\n\
         HMM          A        C\n\
         \x20           m->m     m->i     m->d     i->m     i->i     d->m     d->d\n\
         \x20      1.00000  1.00000\n\
         \x20      0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000\n"
    );

    for node in 1..=leng {
        out.push_str(&format!("{node:>7}  1.00000  1.00000\n"));
        out.push_str("       1.00000  1.00000\n");
        out.push_str("       0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000\n");
    }

    out.push_str("//\n");
    out
}
