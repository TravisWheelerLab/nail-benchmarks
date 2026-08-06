pub mod a {
    const ROUNDS: usize = 5;

    pub trait RoundFunction {
        fn new(seed: u64, rounds: usize) -> Self;
        fn apply(&self, round: usize, x: u64) -> u64;
    }

    #[derive(Clone, Debug)]
    pub struct IntegerMixer {
        round_keys: Vec<u64>,
    }

    impl IntegerMixer {}

    impl RoundFunction for IntegerMixer {
        fn new(seed: u64, rounds: usize) -> Self {
            let mut state = seed;
            let round_keys = (0..rounds).map(|_| splitmix64_next(&mut state)).collect();

            Self { round_keys }
        }

        #[inline]
        fn apply(&self, round: usize, x: u64) -> u64 {
            mix64(x ^ self.round_keys[round])
        }
    }

    #[inline]
    fn mix64(mut x: u64) -> u64 {
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
        x ^= x >> 31;
        x
    }

    #[inline]
    fn splitmix64_next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        mix64(*state)
    }

    pub struct Feistel<R>
    where
        R: RoundFunction,
    {
        round_fn: R,
        bits: usize,
    }

    impl<R> Feistel<R>
    where
        R: RoundFunction,
    {
        pub fn new(bits: usize, seed: u64) -> Feistel<R> {
            assert!(bits.is_multiple_of(2));
            assert!(bits <= 64);
            Feistel {
                round_fn: R::new(seed, ROUNDS),
                bits,
            }
        }

        pub fn permute(&self, x: u64) -> u64 {
            let (mut l, mut r) = Self::split(x, self.bits);

            let n = self.bits >> 1;
            let mask = (1u64 << n) - 1;
            for round in 0..ROUNDS {
                l ^= self.round_fn.apply(round, r) & mask;
                (l, r) = (r, l);
            }
            Self::combine(r, l, self.bits)
        }

        fn split(x: u64, bits: usize) -> (u64, u64) {
            let n = bits >> 1;
            let m = (1u64 << n) - 1;
            let hi = x >> n;
            let lo = x & m;
            (hi, lo)
        }

        fn combine(hi: u64, lo: u64, bits: usize) -> u64 {
            let n = bits >> 1;
            (hi << n) | lo
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
            let z = 64 - n.leading_zeros() as usize;
            let bits = z + (z & 1);

            Permutation {
                n,
                feistel: Feistel::new(bits, seed),
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
}

pub mod b {
    use std::hash::{BuildHasher, Hasher};

    pub struct Feistel<B>
    where
        B: BuildHasher,
    {
        hash_builder: B,
        bits: usize,
        keys: Vec<u64>,
    }

    impl<B> Feistel<B>
    where
        B: BuildHasher,
    {
        pub fn new(hash_builder: B, bits: usize, keys: &[u64]) -> Feistel<B> {
            assert!(bits.is_multiple_of(2));
            assert!(bits <= 64);
            Feistel {
                hash_builder,
                bits,
                keys: Vec::from(keys),
            }
        }

        pub fn encrypt(&self, x: u64) -> u64 {
            let (mut l, mut r) = self.split(x);
            for k in self.keys.iter() {
                l ^= self.hash(*k, r);
                (l, r) = (r, l);
            }
            self.combine(r, l)
        }

        pub fn decrypt(&self, x: u64) -> u64 {
            let (mut l, mut r) = self.split(x);
            for k in self.keys.iter().rev() {
                l ^= self.hash(*k, r);
                (l, r) = (r, l);
            }
            self.combine(r, l)
        }

        fn split(&self, x: u64) -> (u64, u64) {
            let n = self.bits >> 1;
            let m = (1u64 << n) - 1;
            let hi = x >> n;
            let lo = x & m;
            (hi, lo)
        }

        fn combine(&self, hi: u64, lo: u64) -> u64 {
            let n = self.bits >> 1;
            (hi << n) | lo
        }

        fn hash(&self, k: u64, x: u64) -> u64 {
            let mut h: <B as BuildHasher>::Hasher = self.hash_builder.build_hasher();
            h.write_u64(k);
            h.write_u64(x);
            let res = h.finish();
            let n = self.bits >> 1;
            let m = (1u64 << n) - 1;
            res & m
        }
    }

    pub struct Permutation<B>
    where
        B: BuildHasher,
    {
        n: u64,
        feistel: Feistel<B>,
    }

    impl<B> Permutation<B>
    where
        B: BuildHasher,
    {
        pub fn new(n: u64, seed: u64, bob: B) -> Permutation<B> {
            let mut keys = Vec::new();
            let mut k = seed;
            for _i in 0..5 {
                k = bob.hash_one(k);
                keys.push(k);
            }

            // Code assumes an even number of bits. Rounding up
            // increases the constant factor in [`get`] but doesn't
            // alter the big-O complexity.
            let z = 64 - n.leading_zeros() as usize;
            let bits = z + (z & 1);

            Permutation {
                n,
                feistel: Feistel::new(bob, bits, &keys),
            }
        }

        pub fn get(&self, x: u64) -> u64 {
            assert!(x < self.n);
            let mut res = self.feistel.encrypt(x);
            while res >= self.n {
                res = self.feistel.encrypt(res);
            }
            res
        }
    }
}
