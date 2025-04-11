#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use whir::crypto::fields::{Field64, Field64_2};

pub mod constraints;
pub mod expr;
pub mod lagrange;
pub mod partial_sums;
pub mod pcs;

pub type BaseField = Field64;
pub type ExtField = Field64_2;
