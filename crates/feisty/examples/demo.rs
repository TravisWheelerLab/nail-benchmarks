use feisty::{IntegerMixer, Permutation};

type Perm = Permutation<IntegerMixer>;

fn main() {
    let seed = 67779;

    // tiny example
    let perm = Perm::new(16, seed);
    let shuffled: Vec<u64> = (0..16).map(|i| perm.get(i)).collect();
    println!("n=16: {shuffled:?}");

    let mut seen = shuffled.clone();
    seen.sort_unstable();
    assert!(seen.into_iter().eq(0..16), "not a permutation");

    // different seed = different permuation
    let other = Perm::new(16, seed + 1);
    let shuffled: Vec<u64> = (0..16).map(|i| other.get(i)).collect();
    println!("n=16: {shuffled:?}  (seed+1)");

    // the actual use case: large permutations are fast and don't require large memory
    let n = 1_000_000_000;
    let perm = Perm::new(n, seed);
    println!("\nn={n}");
    for i in [0, 1, 2, n / 2, n - 1] {
        println!("  {i:>10} -> {:>10}", perm.get(i));
    }
}
