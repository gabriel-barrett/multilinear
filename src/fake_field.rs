use std::{
    iter::Product,
    ops::{Div, DivAssign, Mul, MulAssign},
};

// This is a fake field, only for benchmarking purposes
use derive_more::{Add, AddAssign, From, Into, Neg, Sub, SubAssign, Sum};

#[derive(
    Add,
    AddAssign,
    Clone,
    Copy,
    Debug,
    Default,
    Eq,
    From,
    Hash,
    Into,
    Neg,
    Ord,
    PartialEq,
    PartialOrd,
    Sub,
    SubAssign,
    Sum,
)]
#[repr(transparent)]
pub struct Fake64(i64);

impl Product for Fake64 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Fake64(iter.map(|x| x.0).product())
    }
}

impl Mul for Fake64 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Fake64(self.0 * rhs.0)
    }
}

impl MulAssign for Fake64 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = Fake64(self.0 * rhs.0);
    }
}

impl Div for Fake64 {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        Fake64(self.0 / rhs.0)
    }
}

impl DivAssign for Fake64 {
    fn div_assign(&mut self, rhs: Self) {
        *self = Fake64(self.0 / rhs.0);
    }
}

impl From<u64> for Fake64 {
    fn from(x: u64) -> Self {
        Fake64(x as i64)
    }
}

impl From<i32> for Fake64 {
    fn from(x: i32) -> Self {
        Fake64(x as i64)
    }
}

impl From<u32> for Fake64 {
    fn from(x: u32) -> Self {
        Fake64(x as i64)
    }
}
