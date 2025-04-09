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

#[cfg(test)]
mod tests {
    use crate::BaseField;

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
