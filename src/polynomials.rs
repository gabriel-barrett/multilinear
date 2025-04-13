use std::ops::Add;

use ark_ff::Field;
use whir::poly_utils::coeffs::CoefficientList;

#[derive(Default, Debug, Clone)]
pub struct LagrangePolynomial<F>(Vec<F>);

impl<F: Field> LagrangePolynomial<F> {
    pub fn new(evals: Vec<F>) -> Self {
        let len = evals.len();
        assert!(len.is_power_of_two());
        Self(evals)
    }

    pub fn to_coefficient_list(self) -> CoefficientList<F> {
        let coeffs = mobius_inversion(self.0);
        CoefficientList::new(coeffs)
    }
}

fn mobius_inversion<F: Field>(evals: Vec<F>) -> Vec<F> {
    let n = evals.len().trailing_zeros(); // Assumes evals.len() == 2^n
    let mut coeffs = evals;

    for i in 0..n {
        for mask in 0..(1 << n) {
            if (mask & (1 << i)) != 0 {
                let tmp = coeffs[mask ^ (1 << i)];
                coeffs[mask] -= tmp;
            }
        }
    }

    coeffs
}

#[derive(Clone, Copy, Debug)]
pub struct SquarePolynomialEval<F> {
    pub one: F,
    pub zero: F,
    pub minus_one: F,
}

impl<F: Field> Add for SquarePolynomialEval<F> {
    type Output = SquarePolynomialEval<F>;
    fn add(self, other: SquarePolynomialEval<F>) -> SquarePolynomialEval<F> {
        SquarePolynomialEval {
            one: self.one + other.one,
            zero: self.zero + other.zero,
            minus_one: self.minus_one + other.minus_one,
        }
    }
}

impl<F: Field> SquarePolynomialEval<F> {
    pub fn to_polynomial(self) -> SquarePolynomial<F> {
        let a = self.zero;
        let b = (self.one - self.minus_one) / F::from(2);
        let c = (self.one + self.minus_one) / F::from(2) - a;
        SquarePolynomial { coeffs: [a, b, c] }
    }

    pub fn evaluate(&self, x: F) -> F {
        self.to_polynomial().evaluate(x)
    }
}

#[derive(Clone, Copy)]
pub struct SquarePolynomial<F> {
    pub coeffs: [F; 3],
}

impl<F: Field> std::fmt::Debug for SquarePolynomial<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [a, b, c] = self.coeffs;
        write!(f, "{} + {}*X + {}*X^2", a, b, c)?;
        Ok(())
    }
}

impl<F: Field> SquarePolynomial<F> {
    pub fn evaluate(&self, x: F) -> F {
        let [a, b, c] = self.coeffs;
        a + b * x + c * x * x
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::BaseField;

    use super::*;
    use whir::poly_utils::multilinear::MultilinearPoint;

    #[test]
    fn translation_test() {
        let f = BaseField::from;
        let coeffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16].map(f);

        let polynomial = CoefficientList::new(coeffs.to_vec());
        let evals = [
            [f(0), f(0), f(0), f(0)],
            [f(0), f(0), f(0), f(1)],
            [f(0), f(0), f(1), f(0)],
            [f(0), f(0), f(1), f(1)],
            [f(0), f(1), f(0), f(0)],
            [f(0), f(1), f(0), f(1)],
            [f(0), f(1), f(1), f(0)],
            [f(0), f(1), f(1), f(1)],
            [f(1), f(0), f(0), f(0)],
            [f(1), f(0), f(0), f(1)],
            [f(1), f(0), f(1), f(0)],
            [f(1), f(0), f(1), f(1)],
            [f(1), f(1), f(0), f(0)],
            [f(1), f(1), f(0), f(1)],
            [f(1), f(1), f(1), f(0)],
            [f(1), f(1), f(1), f(1)],
        ]
        .map(|arr| polynomial.evaluate(&MultilinearPoint(arr.to_vec())))
        .to_vec();
        assert_eq!(coeffs.to_vec(), mobius_inversion(evals));
    }
}
