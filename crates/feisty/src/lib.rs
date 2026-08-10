mod splitmix64;

pub const DEFAULT_ROUNDS: usize = 5;

pub trait RoundFunction {
    fn new(seed: u64, rounds: usize) -> Self;
    fn apply(&self, round: usize, x: u64) -> u64;
}

#[derive(Clone, Debug)]
pub struct IntegerMixer {
    round_keys: Vec<u64>,
}

impl RoundFunction for IntegerMixer {
    fn new(seed: u64, rounds: usize) -> Self {
        Self {
            round_keys: splitmix64::generate(rounds, seed),
        }
    }

    #[inline]
    fn apply(&self, round: usize, x: u64) -> u64 {
        splitmix64::mix(x ^ self.round_keys[round])
    }
}

pub struct Feistel<R>
where
    R: RoundFunction,
{
    round_fn: R,
    half: u32,
    mask: u64,
    rounds: usize,
}

impl<R> Feistel<R>
where
    R: RoundFunction,
{
    pub fn new(bits: usize, seed: u64) -> Feistel<R> {
        Self::with_rounds(bits, seed, DEFAULT_ROUNDS)
    }

    pub fn with_rounds(bits: usize, seed: u64, rounds: usize) -> Feistel<R> {
        assert!(bits.is_multiple_of(2));
        assert!(bits <= 64);

        let half = (bits >> 1) as u32;
        Feistel {
            round_fn: R::new(seed, rounds),
            half,
            mask: (1u64 << half) - 1,
            rounds,
        }
    }

    pub fn permute(&self, x: u64) -> u64 {
        let (mut l, mut r) = self.split(x);

        for round in 0..self.rounds {
            l ^= self.round_fn.apply(round, r) & self.mask;
            (l, r) = (r, l);
        }

        self.combine(r, l)
    }

    #[inline]
    fn split(&self, x: u64) -> (u64, u64) {
        (x >> self.half, x & self.mask)
    }

    #[inline]
    fn combine(&self, hi: u64, lo: u64) -> u64 {
        (hi << self.half) | lo
    }
}

pub struct Permutation<R>
where
    R: RoundFunction,
{
    n: u64,
    feistel: Feistel<R>,
}

impl<R> Permutation<R>
where
    R: RoundFunction,
{
    pub fn new(n: u64, seed: u64) -> Permutation<R> {
        Self::with_rounds(n, seed, DEFAULT_ROUNDS)
    }

    pub fn with_rounds(n: u64, seed: u64, rounds: usize) -> Permutation<R> {
        assert!(n > 0);

        // what: compute the smallest even number
        //       of bits required to store n - 1
        //
        // why:  since we are permuting 0..n, we
        //       never actually need to represent n
        //
        //       we want to use fewer bits since it
        //       constrains the output range of the
        //       Feistel network, which means we need
        //       fewer permute() cycles on average to
        //       produce an output in the correct range
        let z = 64 - (n - 1).leading_zeros() as usize;
        let bits = z.next_multiple_of(2);

        Permutation {
            n,
            feistel: Feistel::with_rounds(bits, seed, rounds),
        }
    }

    pub fn get(&self, x: u64) -> u64 {
        assert!(x < self.n);

        let mut res = self.feistel.permute(x);

        while res >= self.n {
            res = self.feistel.permute(res);
        }

        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_bijection(n: u64, seed: u64) {
        let perm = Permutation::<IntegerMixer>::new(n, seed);
        let mut seen = vec![false; n as usize];

        for x in 0..n {
            let y = perm.get(x);
            assert!(y < n, "n={n} seed={seed}: get({x}) = {y} out of range");
            assert!(!seen[y as usize], "n={n} seed={seed}: {y} seen twice");
            seen[y as usize] = true;
        }

        for x in 0..n {
            assert!(seen[x as usize], "n={n} seed={seed}: {x} missed");
        }
    }

    #[test]
    fn permutes_every_small_n() {
        for n in 1..=1024 {
            for seed in [0, 1, 67779] {
                assert_bijection(n, seed);
            }
        }
    }

    #[test]
    fn permutes_powers_of_two() {
        for k in 0..=16 {
            assert_bijection(1 << k, 67779);
        }
    }

    #[test]
    fn single_element_is_identity() {
        assert_eq!(Permutation::<IntegerMixer>::new(1, 67779).get(0), 0);
    }

    #[test]
    fn seeds_give_different_permutations() {
        let n = 1000;
        let a = Permutation::<IntegerMixer>::new(n, 1);
        let b = Permutation::<IntegerMixer>::new(n, 2);
        assert!((0..n).any(|x| a.get(x) != b.get(x)));
    }

    #[test]
    fn is_deterministic() {
        let n = 1000;
        let a = Permutation::<IntegerMixer>::new(n, 67779);
        let b = Permutation::<IntegerMixer>::new(n, 67779);
        assert!((0..n).all(|x| a.get(x) == b.get(x)));
    }

    #[test]
    #[should_panic]
    fn rejects_empty_range() {
        Permutation::<IntegerMixer>::new(0, 67779);
    }

    #[test]
    #[should_panic]
    fn rejects_out_of_range_input() {
        Permutation::<IntegerMixer>::new(10, 67779).get(10);
    }
}
