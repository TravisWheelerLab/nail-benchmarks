use crate::splitmix64;

pub const DEFAULT_ROUNDS: usize = 5;

pub trait RoundFunction {
    fn new(seed: u64, rounds: usize) -> Self;
    fn apply(&self, round: usize, x: u64) -> u64;
}

#[derive(Clone, Debug)]
pub struct IntegerMixer {
    keys: [u64; IntegerMixer::MAX_ROUNDS],
}

impl IntegerMixer {
    pub const MAX_ROUNDS: usize = 16;
}

impl RoundFunction for IntegerMixer {
    fn new(seed: u64, rounds: usize) -> Self {
        assert!(rounds <= Self::MAX_ROUNDS);

        let mut keys = [0u64; Self::MAX_ROUNDS];
        let mut state = seed;
        for key in keys.iter_mut().take(rounds) {
            *key = splitmix64::next(&mut state);
        }

        Self { keys }
    }

    #[inline]
    fn apply(&self, round: usize, x: u64) -> u64 {
        // new() caps rounds, so the clamp never changes the index.
        // it's here to tell the compiler there is nothing to bounds check.
        splitmix64::mix(x ^ self.keys[round.min(Self::MAX_ROUNDS - 1)])
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
