// references:
//
// - SplitMix64: the algorithm implemented here, a single-number state driving
//   a bit-mixing process.
//   Guy Steele, Doug Lea and Christine Flood, "Fast Splittable Pseudorandom
//   Number Generators", OOPSLA 2014
//
// - MurmurHash3: where the mixer's shape came from.
//   Austin Appleby, 2011
//
// - Mix13: the constants and shift amounts in [`mix`], found by search.
//   David Stafford, "Better Bit Mixing", 2011
//
// - splitmix64.c: the reference implementation [`next`] is tested against.
//   Sebastiano Vigna, 2015

/// SplitMix64's golden-ratio increment.
const GOLDEN_GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;

/// Stafford's Mix13 multipliers. The shift amounts in [`mix`] came out of the
/// same search, so do not change one without the others.
const MIX_MUL_A: u64 = 0xbf58_476d_1ce4_e5b9;
const MIX_MUL_B: u64 = 0x94d0_49bb_1331_11eb;

/// Mixes the bits of `x` into a completely different number.
#[inline]
pub fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(MIX_MUL_A);
    x ^= x >> 27;
    x = x.wrapping_mul(MIX_MUL_B);
    x ^= x >> 31;
    x
}

/// Steps the state forward and hands back a mixed copy of it.
#[inline]
pub fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(GOLDEN_GAMMA);
    mix(*state)
}

/// Generate `n` random numbers using `seed`.
///
/// This does NOT use any global state, so repeated calls with
/// the same arguments will always return the same numbers.
pub fn generate(n: usize, seed: u64) -> Vec<u64> {
    let mut state = seed;
    (0..n).map(|_| next(&mut state)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first five numbers Vigna's splitmix64.c produces starting from 0.
    #[test]
    fn matches_reference_vectors() {
        let correct = [
            0xE220_A839_7B1D_CDAFu64,
            0x6E78_9E6A_A1B9_65F4,
            0x06C4_5D18_8009_454F,
            0xF88B_B8A8_724C_81EC,
            0x1B39_896A_51A8_749B,
        ];

        let mut state = 0u64;
        for (i, &w) in correct.iter().enumerate() {
            assert_eq!(next(&mut state), w, "output {i}");
        }
    }

    #[test]
    fn mix_of_zero_is_zero() {
        assert_eq!(mix(0), 0);
    }
}
