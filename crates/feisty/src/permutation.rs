use crate::feistel::{Feistel, IntegerMixer, RoundFunction, DEFAULT_ROUNDS};
use crate::Index;

pub struct Permutation<T, R = IntegerMixer>
where
    T: Index,
    R: RoundFunction,
{
    n: T,
    feistel: Feistel<R>,
}

impl<T, R> Permutation<T, R>
where
    T: Index,
    R: RoundFunction,
{
    pub fn with_rounds(n: T, seed: u64, rounds: usize) -> Permutation<T, R> {
        let size = n.to_u64();
        assert!(size > 0);

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
        let z = 64 - (size - 1).leading_zeros() as usize;
        let bits = z.next_multiple_of(2);

        Permutation {
            n,
            feistel: Feistel::with_rounds(bits, seed, rounds),
        }
    }

    pub fn get(&self, x: T) -> T {
        let n = self.n.to_u64();
        let x = x.to_u64();
        assert!(x < n);

        let mut res = self.feistel.permute(x);

        while res >= n {
            res = self.feistel.permute(res);
        }

        T::from_u64(res)
    }
}

// new() lives here on its own so the mixer never has to be named.
// nothing about R shows up in the arguments, so anything spelled
// Permutation::with_rounds has to say which mixer it wants.
impl<T> Permutation<T, IntegerMixer>
where
    T: Index,
{
    pub fn new(n: T, seed: u64) -> Permutation<T, IntegerMixer> {
        Self::with_rounds(n, seed, DEFAULT_ROUNDS)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The loop runs in u64 so one helper covers every index type.
    fn assert_bijection<T: Index>(n: u64, seed: u64) {
        let perm = Permutation::<T>::new(T::from_u64(n), seed);
        let mut seen = vec![false; n as usize];

        for x in 0..n {
            let y = perm.get(T::from_u64(x)).to_u64();
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
                assert_bijection::<u64>(n, seed);
            }
        }
    }

    #[test]
    fn permutes_powers_of_two() {
        for k in 0..=16 {
            assert_bijection::<u64>(1 << k, 67779);
        }
    }

    #[test]
    fn permutes_for_every_index_type() {
        assert_bijection::<u8>(255, 67779);
        assert_bijection::<u16>(1000, 67779);
        assert_bijection::<u32>(1000, 67779);
        assert_bijection::<u64>(1000, 67779);
        assert_bijection::<usize>(1000, 67779);
    }

    #[test]
    fn handles_the_largest_domain_of_each_type() {
        Permutation::new(u8::MAX, 67779);
        Permutation::new(u16::MAX, 67779);
        Permutation::new(u32::MAX, 67779);
        Permutation::new(u64::MAX, 67779);
        Permutation::new(usize::MAX, 67779);
    }

    /// Values recorded before the crate went generic; they must not drift.
    #[test]
    fn matches_baseline() {
        let perm = Permutation::new(1_000_000usize, 67779);
        let got: Vec<usize> = (0..5).map(|x| perm.get(x)).collect();
        assert_eq!(got, [693430, 65201, 255456, 22500, 237480]);
    }

    #[test]
    fn single_element_is_identity() {
        assert_eq!(Permutation::new(1u64, 67779).get(0), 0);
    }

    #[test]
    fn seeds_give_different_permutations() {
        let n = 1000u64;
        let a = Permutation::new(n, 1);
        let b = Permutation::new(n, 2);
        assert!((0..n).any(|x| a.get(x) != b.get(x)));
    }

    #[test]
    fn is_deterministic() {
        let n = 1000u64;
        let a = Permutation::new(n, 67779);
        let b = Permutation::new(n, 67779);
        assert!((0..n).all(|x| a.get(x) == b.get(x)));
    }

    #[test]
    #[should_panic]
    fn rejects_empty_range() {
        Permutation::new(0u64, 67779);
    }

    #[test]
    #[should_panic]
    fn rejects_out_of_range_input() {
        Permutation::new(10u64, 67779).get(10);
    }

    #[test]
    #[should_panic]
    fn rejects_too_many_rounds() {
        Permutation::<u64>::with_rounds(1000, 67779, IntegerMixer::MAX_ROUNDS + 1);
    }
}
