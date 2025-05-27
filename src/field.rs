use std::{
    iter::{Product, Sum},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use ark_ff::{Fp128, MontBackend, MontConfig};
use serde::de::Error as DeError;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

pub trait Field = Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Neg<Output = Self>
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
    + Sized
    + Copy
    + Sum
    + Product
    + From<i64>
    + Eq;

#[derive(MontConfig)]
#[modulus = "340282366920938463463374557953744961537"]
#[generator = "3"]
pub struct FrConfig128;
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Field128(pub Fp128<MontBackend<FrConfig128, 2>>);

impl AsRef<[u8]> for Field128 {
    fn as_ref(&self) -> &[u8] {
        let data = &self.0 .0 .0;
        unsafe {
            std::slice::from_raw_parts(
                data.as_ptr() as *const u8,
                data.len() * std::mem::size_of::<u64>(),
            )
        }
    }
}

impl Serialize for Field128 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let bytes = self.as_ref();
        serializer.serialize_bytes(bytes)
    }
}

impl<'de> Deserialize<'de> for Field128 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let bytes: &[u8] = Deserialize::deserialize(deserializer)?;
        if bytes.len() != 16 {
            return Err(D::Error::custom("Invalid byte length for Field128"));
        }
        let mut array = [0u8; 16];
        array.copy_from_slice(bytes);
        let value = u128::from_le_bytes(array);
        Ok(Field128(Fp128::from(value)))
    }
}

// Implement Add
impl Add for Field128 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Field128(self.0 + rhs.0)
    }
}

impl AddAssign for Field128 {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl Sub for Field128 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Field128(self.0 - rhs.0)
    }
}

impl SubAssign for Field128 {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl Neg for Field128 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Field128(-self.0)
    }
}

impl Mul for Field128 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Field128(self.0 * rhs.0)
    }
}

impl MulAssign for Field128 {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}

impl Div for Field128 {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        Field128(self.0 / rhs.0)
    }
}

impl DivAssign for Field128 {
    fn div_assign(&mut self, rhs: Self) {
        self.0 /= rhs.0;
    }
}

impl Sum for Field128 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Field128(Fp128::from(0)), |a, b| a + b)
    }
}

impl Product for Field128 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Field128(Fp128::from(1)), |a, b| a * b)
    }
}

impl From<u128> for Field128 {
    fn from(val: u128) -> Self {
        Field128(Fp128::from(val))
    }
}

impl From<i32> for Field128 {
    fn from(val: i32) -> Self {
        Field128(Fp128::from(val))
    }
}

impl From<i64> for Field128 {
    fn from(val: i64) -> Self {
        Field128(Fp128::from(val))
    }
}

impl std::fmt::Display for Field128 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
