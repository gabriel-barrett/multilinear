use whir::crypto::fields::{Field64, Field64_2};

use rand::{rngs::ThreadRng, Rng};

pub fn random_base(rng: &mut ThreadRng) -> BaseField {
    BaseField::from(rng.random::<u64>())
}

pub type BaseField = Field64;
pub type ExtField = Field64_2;
