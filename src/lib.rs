#![feature(trait_alias)]
pub mod constraint_system;
pub mod field;
pub mod fri;
pub mod merkle_tree;
pub mod ntt;
pub mod polynomials;
pub mod transcript;

#[cfg(test)]
pub mod p3;

#[macro_export]
macro_rules! benchmark {
    ($msg:expr, $expr:expr) => {{
        let now = std::time::Instant::now();
        let result = std::hint::black_box($expr);
        println!("{}{:?}", $msg, now.elapsed());
        result
    }};
}
