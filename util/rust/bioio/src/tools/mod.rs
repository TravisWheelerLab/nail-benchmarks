pub mod hmmer;
pub mod mmseqs;
pub mod nail;

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
