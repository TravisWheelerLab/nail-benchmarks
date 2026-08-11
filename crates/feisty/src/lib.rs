mod feistel;
mod permutation;
mod splitmix64;

pub use feistel::{DEFAULT_ROUNDS, Feistel, IntegerMixer, RoundFunction};
pub use permutation::Permutation;

// Index::to_u64 casts, so a usize wider than u64 would lose the top bits
const _: () = assert!(usize::BITS <= 64);

mod sealed {
    pub trait Sealed {}

    impl Sealed for u8 {}
    impl Sealed for u16 {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
    impl Sealed for usize {}
}

/// An unsigned integer a [`Permutation`] can take in and hand back.
pub trait Index: Copy + sealed::Sealed {
    fn to_u64(self) -> u64;
    fn from_u64(x: u64) -> Self;
}

macro_rules! impl_index {
    ($($t:ty)*) => {$(
        impl Index for $t {
            #[inline]
            fn to_u64(self) -> u64 {
                self as u64
            }

            #[inline]
            fn from_u64(x: u64) -> Self {
                x as Self
            }
        }
    )*};
}

impl_index!(u8 u16 u32 u64 usize);
