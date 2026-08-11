use feisty::Permutation;

fn main() {
    let seed = 67779;

    // tiny example
    let perm = Permutation::new(16usize, seed);
    let shuffled: Vec<usize> = (0..16).map(|i| perm.get(i)).collect();
    println!("n=16: {shuffled:?}");

    let mut seen = shuffled.clone();
    seen.sort_unstable();
    assert!(seen.into_iter().eq(0..16), "not a permutation");

    // different seed = different permuation
    let other = Permutation::new(16usize, seed + 1);
    let shuffled: Vec<usize> = (0..16).map(|i| other.get(i)).collect();
    println!("n=16: {shuffled:?}  (seed+1)");

    // whatever integer type you hand it is the one you get back
    let small = Permutation::new(200u8, seed);
    let wide = Permutation::new(1u64 << 40, seed);
    println!("\nu8:  {} -> {}", 7u8, small.get(7));
    println!("u64: {} -> {}", 7u64, wide.get(7));

    // the actual use case: large permutations are fast and don't require large memory
    let n = 1_000_000_000usize;
    let perm = Permutation::new(n, seed);
    println!("\nn={n}");
    for i in [0, 1, 2, n / 2, n - 1] {
        println!("  {i:>10} -> {:>10}", perm.get(i));
    }
}
