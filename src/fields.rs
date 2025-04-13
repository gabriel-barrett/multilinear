//! An implementation of a 64-bit STARK-friendly prime field with modulus `2^64 -
//! 2^32 + 1`. The implementation follows <https://eprint.iacr.org/2022/274.pdf>
//! and the code for the majority of functions was taken and adapted from
//! <https://github.com/novifinancial/winterfell>
//!
//! This field and its implementation has many attractive properties:
//! * Multiplication of two 32-bit values does not overflow field modulus.
//! * Field arithmetic in this field can be implemented using a few 32-bit
//!   addition, subtractions, and shifts.
//! * 8 is the 64th root of unity which opens up potential for optimized FFT
//!   implementations.

use ark_ff::{fields::Fp64, BigInt, Fp, Fp2, Fp2Config, PrimeField, SqrtPrecomputation, Zero};
use core::marker::PhantomData;
use rand::{rngs::ThreadRng, Rng};

pub fn random_base(rng: &mut ThreadRng) -> BaseField {
    BaseField::from(rng.random::<u64>())
}

/// An optimized implementation of a 64-bit prime field with modulus `2^64 -
/// 2^32 + 1`
pub type BaseField = Fp64<FpParams>;

/// Extension field
pub type ExtField = Fp2<F2Config64>;

pub struct F2Config64;
impl Fp2Config for F2Config64 {
    type Fp = Fp64<FpParams>;

    const NONRESIDUE: Self::Fp = Fp::<FpParams, 1>(BigInt([7]), PhantomData);

    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        // Fq(7)**(((q^0) - 1) / 2)
        Fp::<FpParams, 1>(BigInt([1]), PhantomData),
        // Fq(7)**(((q^1) - 1) / 2)
        Fp::<FpParams, 1>(BigInt([18446744069414584320]), PhantomData),
    ];
}

/// Field modulus `p = 2^64 - 2^32 + 1`
const MODULUS: u64 = 18_446_744_069_414_584_321;

/// Square of auxiliary modulus R for Montgomery reduction `R2 ≡ (2^64)^2 mod p`
const R2: u64 = 18_446_744_065_119_617_025;

pub struct FpParams;
impl ark_ff::FpConfig<1> for FpParams {
    const MODULUS: ark_ff::BigInt<1> = BigInt([MODULUS]);
    const GENERATOR: Fp64<Self> = into_mont(7);
    const ZERO: Fp64<Self> = into_mont(0);
    const ONE: Fp64<Self> = into_mont(1);
    const TWO_ADICITY: u32 = 32;
    const TWO_ADIC_ROOT_OF_UNITY: Fp64<Self> = into_mont(1_753_635_133_440_165_772);
    const SQRT_PRECOMP: Option<ark_ff::SqrtPrecomputation<Fp64<Self>>> =
        Some(SqrtPrecomputation::TonelliShanks {
            two_adicity: Self::TWO_ADICITY,
            quadratic_nonresidue_to_trace: Self::TWO_ADIC_ROOT_OF_UNITY,
            trace_of_modulus_minus_one_div_two: &<Fp64<Self>>::TRACE_MINUS_ONE_DIV_TWO.0,
        });

    fn add_assign(a: &mut Fp64<Self>, b: &Fp64<Self>) {
        // We compute a + b = a - (p - b).
        let (x1, c1) = (a.0).0[0].overflowing_sub(MODULUS - (b.0).0[0]);
        let adj = 0u32.wrapping_sub(u32::from(c1));
        (a.0).0[0] = x1.wrapping_sub(u64::from(adj));
    }

    fn sub_assign(a: &mut Fp64<Self>, b: &Fp64<Self>) {
        let (x1, c1) = (a.0).0[0].overflowing_sub((b.0).0[0]);
        let adj = 0u32.wrapping_sub(u32::from(c1));
        (a.0).0[0] = x1.wrapping_sub(u64::from(adj));
    }

    fn double_in_place(a: &mut Fp64<Self>) {
        Self::add_assign(a, &a.clone());
    }

    fn mul_assign(a: &mut Fp64<Self>, b: &Fp64<Self>) {
        (a.0).0[0] = mont_red(u128::from((a.0).0[0]) * u128::from((b.0).0[0]));
    }

    fn sum_of_products<const T: usize>(a: &[Fp64<Self>; T], b: &[Fp64<Self>; T]) -> Fp64<Self> {
        a.iter().zip(b).map(|(&a, b)| a * b).sum()
    }

    fn square_in_place(a: &mut Fp64<Self>) {
        let temp = *a;
        Self::mul_assign(a, &temp);
    }

    fn inverse(a: &Fp64<Self>) -> Option<Fp64<Self>> {
        if a.is_zero() {
            None
        } else {
            let a = (a.0).0[0];
            let t2 = exp_acc::<1>(a, a);
            let t3 = exp_acc::<1>(t2, a);
            let t6 = exp_acc::<3>(t3, t3);
            let t12 = exp_acc::<6>(t6, t6);
            let t24 = exp_acc::<12>(t12, t12);
            let t30 = exp_acc::<6>(t24, t6);
            let t31 = exp_acc::<1>(t30, a);
            let t63 = exp_acc::<32>(t31, t31);
            let inv = exp_acc::<1>(t63, a);
            Some(ark_ff::Fp(BigInt([inv]), PhantomData))
        }
    }

    fn from_bigint(other: ark_ff::BigInt<1>) -> Option<Fp64<Self>> {
        let inner = other.0[0];
        if inner.is_zero() {
            Some(Self::ZERO)
        } else if inner < MODULUS {
            Some(into_mont(other.0[0]))
        } else {
            None
        }
    }

    fn into_bigint(other: Fp64<Self>) -> ark_ff::BigInt<1> {
        BigInt([mont_red(u128::from((other.0).0[0]))])
    }

    fn neg_in_place(a: &mut Fp64<Self>) {
        let mut tmp = Self::ZERO;
        Self::sub_assign(&mut tmp, a);
        a.0 = tmp.0;
    }
}

/// Converts a value into Montgomery representation
#[inline]
const fn into_mont(value: u64) -> Fp64<FpParams> {
    ark_ff::Fp(BigInt([mont_red(value as u128 * R2 as u128)]), PhantomData)
}

/// Performs Montgomery reduction
#[inline]
const fn mont_red(x: u128) -> u64 {
    // See reference above for a description of the following implementation.
    #[allow(clippy::cast_possible_truncation)]
    let xl = x as u64;
    let xh = (x >> 64) as u64;
    let (a, overflow) = xl.overflowing_add(xl << 32);
    let b = a.wrapping_sub(a >> 32).wrapping_sub(overflow as u64);
    let (r, underflow) = xh.overflowing_sub(b);
    r.wrapping_sub(0u32.wrapping_sub(underflow as u32) as u64)
}

/// Squares `base` N times and multiplies the result by the tail value.
#[inline]
const fn exp_acc<const N: usize>(base: u64, tail: u64) -> u64 {
    let mut result = base;
    let mut i = 0;
    while i < N {
        result = mont_red(result as u128 * result as u128);
        i += 1;
    }
    mont_red(result as u128 * tail as u128)
}
