use anyhow::{ensure, Context, Result};
use ark_ff::AdditiveGroup;

use crate::{
    constraints::BQCS,
    polynomials::{SquarePolynomial, SquarePolynomialEval},
    BaseField,
};
use rand::{rngs::ThreadRng, Rng};

impl<const I: usize, const J: usize, const K: usize> BQCS<I, J, K>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
    [(); 1 << K]:,
{
    // Okay as it is linear in the height, but it is still width squared,
    // even if the constraint matrix usually has a sparse representation.
    // Needs further optimizations.
    pub fn full_sum(&self) -> BaseField {
        let mut acc = BaseField::ZERO;
        for i in 0..1 << I {
            for j1 in 0..1 << J {
                for j2 in 0..1 << J {
                    acc += self.evaluate_hypercube(&[], &[], i, &[], j1, &[], j2);
                }
            }
        }
        acc
    }

    pub fn partial_sum_rows(&self, row1: &[BaseField], row2: &[BaseField]) -> BaseField {
        let len = row1.len();
        assert_eq!(row2.len(), len);
        assert!(len > 0);
        assert!(len <= I);
        let mut acc = BaseField::ZERO;
        for i in 0..1 << (I - len) {
            for j1 in 0..1 << J {
                for j2 in 0..1 << J {
                    acc += self.evaluate_hypercube_col(row1, row2, i, j1, j2);
                }
            }
        }
        acc
    }

    pub fn partial_sum_cols(&self, col1: &[BaseField], col2: &[BaseField]) -> BaseField {
        let len = col1.len();
        assert_eq!(col2.len(), len);
        assert!(len > 0);
        assert!(len <= I);
        let mut acc = BaseField::ZERO;
        for i in 0..1 << I {
            for j1 in 0..1 << (J - len) {
                for j2 in 0..1 << (J - len) {
                    acc += self.evaluate_hypercube_row(i, col1, j1, col2, j2);
                }
            }
        }
        acc
    }

    pub fn partial_sum_cols_fixed_row(
        &self,
        row1: &[BaseField; I],
        row2: &[BaseField; I],
        col1: &[BaseField],
        col2: &[BaseField],
    ) -> BaseField {
        let len = col1.len();
        assert_eq!(col2.len(), len);
        assert!(len > 0);
        assert!(len <= J);
        let mut acc = BaseField::ZERO;
        for j1 in 0..1 << (J - len) {
            for j2 in 0..1 << (J - len) {
                acc += self.evaluate_hypercube(row1, row2, 0, col1, j1, col2, j2);
            }
        }
        acc
    }

    pub fn partial_sum_rows_fixed_col(
        &self,
        row1: &[BaseField],
        row2: &[BaseField],
        col1: &[BaseField; J],
        col2: &[BaseField; J],
    ) -> BaseField {
        let len = row1.len();
        assert_eq!(row2.len(), len);
        assert!(len > 0);
        assert!(len <= I);
        let mut acc = BaseField::ZERO;
        for i in 0..1 << (I - len) {
            acc += self.evaluate_hypercube(row1, row2, i, col1, 0, col2, 0);
        }
        acc
    }

    pub fn generate_partial_polynomials(&self) -> Vec<(BaseField, SquarePolynomial<BaseField>)> {
        let mut rng = rand::rng();
        let mut acc = Vec::with_capacity(2 * I + 2 * J);

        let mut index_j1 = Vec::with_capacity(J);
        let mut index_j2 = Vec::with_capacity(J);
        for _ in 0..J {
            self.push_polynomials(
                &mut index_j1,
                &mut index_j2,
                &mut rng,
                &mut acc,
                |index1, index2| self.partial_sum_cols(index1, index2),
            );
        }

        let mut iter = acc.iter().map(|(r, _)| r);
        let cols = [(); J].map(|_| (*iter.next().unwrap(), *iter.next().unwrap()));
        let col1 = &cols.map(|(col, _)| col);
        let col2 = &cols.map(|(_, col)| col);
        let mut index_i1 = Vec::with_capacity(I);
        let mut index_i2 = Vec::with_capacity(I);
        for _ in 0..I {
            self.push_polynomials(
                &mut index_i1,
                &mut index_i2,
                &mut rng,
                &mut acc,
                |index1, index2| self.partial_sum_rows_fixed_col(index1, index2, col1, col2),
            );
        }

        acc
    }

    pub fn generate_partial_polynomials_transposed(
        &self,
    ) -> Vec<(BaseField, SquarePolynomial<BaseField>)> {
        let mut rng = rand::rng();
        let mut acc = Vec::with_capacity(2 * I + 2 * J);

        let mut index_i1 = Vec::with_capacity(I);
        let mut index_i2 = Vec::with_capacity(I);
        for _ in 0..I {
            self.push_polynomials(
                &mut index_i1,
                &mut index_i2,
                &mut rng,
                &mut acc,
                |index1, index2| self.partial_sum_rows(index1, index2),
            );
        }

        let mut iter = acc.iter().map(|(r, _)| r);
        let rows = [(); I].map(|_| (*iter.next().unwrap(), *iter.next().unwrap()));
        let row1 = &rows.map(|(row, _)| row);
        let row2 = &rows.map(|(_, row)| row);
        let mut index_j1 = Vec::with_capacity(J);
        let mut index_j2 = Vec::with_capacity(J);
        for _ in 0..J {
            self.push_polynomials(
                &mut index_j1,
                &mut index_j2,
                &mut rng,
                &mut acc,
                |index1, index2| self.partial_sum_cols_fixed_row(row1, row2, index1, index2),
            );
        }

        acc
    }

    fn push_polynomials(
        &self,
        index1: &mut Vec<BaseField>,
        index2: &mut Vec<BaseField>,
        rng: &mut ThreadRng,
        acc: &mut Vec<(BaseField, SquarePolynomial<BaseField>)>,
        partial_sum: impl Fn(&[BaseField], &[BaseField]) -> BaseField,
    ) {
        let f = BaseField::from;
        index1.push(f(-1));
        index2.push(f(0));
        let minus_one0 = partial_sum(index1, index2);
        let minus_one1 = partial_sum(index1, modify_last(index2, f(1)));
        let zero1 = partial_sum(modify_last(index1, f(0)), index2);
        let zero0 = partial_sum(index1, modify_last(index2, f(0)));
        let one0 = partial_sum(modify_last(index1, f(1)), index2);
        let one1 = partial_sum(index1, modify_last(index2, f(1)));
        let pol10 = SquarePolynomialEval {
            minus_one: minus_one0,
            zero: zero0,
            one: one0,
        };
        let pol11 = SquarePolynomialEval {
            minus_one: minus_one1,
            zero: zero1,
            one: one1,
        };
        let pol1 = pol10 + pol11;
        let r1 = BaseField::from(rng.random::<u64>());
        let pol2 = SquarePolynomialEval {
            minus_one: partial_sum(modify_last(index1, r1), modify_last(index2, f(-1))),
            zero: pol10.evaluate(r1),
            one: pol11.evaluate(r1),
        };
        let r2 = BaseField::from(rng.random::<u64>());
        modify_last(index2, r2);
        acc.push((r1, pol1.to_polynomial()));
        acc.push((r2, pol2.to_polynomial()));
    }

    #[allow(dead_code)]
    fn push_polynomials_no_opt(
        &self,
        index1: &mut Vec<BaseField>,
        index2: &mut Vec<BaseField>,
        rng: &mut ThreadRng,
        acc: &mut Vec<(BaseField, SquarePolynomial<BaseField>)>,
        partial_sum: impl Fn(&[BaseField], &[BaseField]) -> BaseField,
    ) {
        let f = BaseField::from;
        index1.push(f(-1));
        index2.push(f(0));
        let minus_one =
            partial_sum(index1, index2) + partial_sum(index1, modify_last(index2, f(1)));
        let zero = partial_sum(modify_last(index1, f(0)), index2)
            + partial_sum(index1, modify_last(index2, f(0)));
        let one = partial_sum(modify_last(index1, f(1)), index2)
            + partial_sum(index1, modify_last(index2, f(1)));
        let pol1 = SquarePolynomialEval {
            minus_one,
            zero,
            one,
        };
        let r1 = BaseField::from(rng.random::<u64>());
        let pol2 = SquarePolynomialEval {
            minus_one: partial_sum(modify_last(index1, r1), modify_last(index2, f(-1))),
            zero: partial_sum(index1, modify_last(index2, f(0))),
            one: partial_sum(index1, modify_last(index2, f(1))),
        };
        let r2 = BaseField::from(rng.random::<u64>());
        modify_last(index2, r2);
        acc.push((r1, pol1.to_polynomial()));
        acc.push((r2, pol2.to_polynomial()));
    }

    pub fn verify_sumcheck(
        &self,
        pols: &[(BaseField, SquarePolynomial<BaseField>)],
        sum: BaseField,
    ) -> Result<()> {
        let mut r_iter = pols.iter().map(|(r, _)| r);
        let cols = [(); J].map(|_| (*r_iter.next().unwrap(), *r_iter.next().unwrap()));
        let rows = [(); I].map(|_| (*r_iter.next().unwrap(), *r_iter.next().unwrap()));
        let row1 = rows.map(|(row, _)| row);
        let row2 = rows.map(|(_, row)| row);
        let col1 = cols.map(|(col, _)| col);
        let col2 = cols.map(|(_, col)| col);
        self.verify_sumcheck_constraints(pols, sum, &row1, &row2, &col1, &col2)
    }

    pub fn verify_sumcheck_transposed(
        &self,
        pols: &[(BaseField, SquarePolynomial<BaseField>)],
        sum: BaseField,
    ) -> Result<()> {
        let mut r_iter = pols.iter().map(|(r, _)| r);
        let rows = [(); I].map(|_| (*r_iter.next().unwrap(), *r_iter.next().unwrap()));
        let cols = [(); J].map(|_| (*r_iter.next().unwrap(), *r_iter.next().unwrap()));
        let row1 = rows.map(|(row, _)| row);
        let row2 = rows.map(|(_, row)| row);
        let col1 = cols.map(|(col, _)| col);
        let col2 = cols.map(|(_, col)| col);
        self.verify_sumcheck_constraints(pols, sum, &row1, &row2, &col1, &col2)
    }

    fn verify_sumcheck_constraints(
        &self,
        pols: &[(BaseField, SquarePolynomial<BaseField>)],
        sum: BaseField,
        row1: &[BaseField; I],
        row2: &[BaseField; I],
        col1: &[BaseField; J],
        col2: &[BaseField; J],
    ) -> Result<()> {
        assert_eq!(pols.len(), 2 * I + 2 * J);
        let f = BaseField::from;
        let mut iter = pols.iter();
        let (mut r, mut pol) = iter.next().context("Expects at least one polynomial")?;
        ensure!(
            sum == pol.evaluate(f(0)) + pol.evaluate(f(1)),
            "Does not sumcheck"
        );
        for (i, (r_next, pol_next)) in iter.enumerate() {
            ensure!(
                pol.evaluate(r) == pol_next.evaluate(f(0)) + pol_next.evaluate(f(1)),
                "Polynomial {i} does not satisfy equation"
            );
            r = *r_next;
            pol = *pol_next;
        }
        ensure!(
            self.evaluate(row1, row2, col1, col2) == pol.evaluate(r),
            "Does not match polynomial evaluation"
        );
        Ok(())
    }
}

// Auxiliary function
fn modify_last(vec: &mut [BaseField], a: BaseField) -> &mut [BaseField] {
    let len = vec.len();
    vec[len - 1] = a;
    vec
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use crate::{
        constraints::{ConstraintSet, Trace},
        expr::Expr,
    };

    const I: usize = 4;
    const J: usize = 2;
    const K: usize = 1;

    fn pythagorean_trace() -> Trace<I, J> {
        let f = BaseField::from;
        // The first three are pythagorean triples
        // The fourth is the sum of the first two
        let rows = Box::new([
            [3, 4, 5, 7].map(f),
            [5, 12, 13, 17].map(f),
            [8, 15, 17, 23].map(f),
            [7, 24, 25, 31].map(f),
            [20, 21, 29, 41].map(f),
            [12, 35, 37, 47].map(f),
            [9, 40, 41, 49].map(f),
            [28, 45, 53, 73].map(f),
            [11, 60, 61, 71].map(f),
            [16, 63, 65, 79].map(f),
            [33, 56, 65, 89].map(f),
            [48, 55, 73, 103].map(f),
            [13, 84, 85, 97].map(f),
            [36, 77, 85, 113].map(f),
            [39, 80, 89, 119].map(f),
            [65, 72, 97, 137].map(f),
        ]);
        Trace::<I, J> { rows }
    }

    fn pythagorean_set() -> ConstraintSet<J, K> {
        let var = Expr::<J>::var;
        let expr1 = var(0) * var(0) + var(1) * var(1) - var(2) * var(2);
        let expr2 = (var(0) + var(1)) * (var(0) + var(1)) - var(3) * var(3);
        let expressions = [expr1, expr2].map(|x| x.to_quadratic().unwrap()).into();
        ConstraintSet::<J, K> { expressions }
    }

    fn pythagorean_cs() -> BQCS<I, J, K> {
        let trace = pythagorean_trace();
        let set = pythagorean_set();
        BQCS::new(trace, set)
    }

    #[test]
    fn full_sum_test() {
        let f = BaseField::from;
        let system = pythagorean_cs();
        assert_eq!(system.full_sum(), f(0));
    }

    #[test]
    fn sumcheck_test() {
        let system = pythagorean_cs();
        let pols = system.generate_partial_polynomials();
        println!("{pols:#?}");
        let f = BaseField::from;
        system.verify_sumcheck(&pols, f(0)).unwrap();
    }

    #[test]
    fn sumcheck_high_bench() {
        let mut rows = pythagorean_trace().rows.to_vec();
        const TOTAL_HEIGHT: usize = 10;
        for _ in 0..TOTAL_HEIGHT - I {
            rows.extend(rows.clone());
        }
        assert_eq!(rows.len(), 1 << TOTAL_HEIGHT);
        let trace = Trace::<TOTAL_HEIGHT, J> {
            rows: rows.try_into().unwrap(),
        };
        let set = pythagorean_set();
        let system = BQCS::new(trace, set);
        println!("GENERATING POLYNOMIAL");
        let now = Instant::now();
        let pols = system.generate_partial_polynomials();
        println!("Generation took {:?}", now.elapsed());
        let f = BaseField::from;
        println!("VERIFYING");
        let now = Instant::now();
        system.verify_sumcheck(&pols, f(0)).unwrap();
        println!("Verification took {:?}", now.elapsed());
    }

    #[test]
    fn sumcheck_high_transposed_bench() {
        let mut rows = pythagorean_trace().rows.to_vec();
        const TOTAL_HEIGHT: usize = 10;
        for _ in 0..TOTAL_HEIGHT - I {
            rows.extend(rows.clone());
        }
        assert_eq!(rows.len(), 1 << TOTAL_HEIGHT);
        let trace = Trace::<TOTAL_HEIGHT, J> {
            rows: rows.try_into().unwrap(),
        };
        let set = pythagorean_set();
        let system = BQCS::new(trace, set);
        println!("GENERATING POLYNOMIAL");
        let now = Instant::now();
        let pols = system.generate_partial_polynomials_transposed();
        println!("Generation took {:?}", now.elapsed());
        let f = BaseField::from;
        println!("VERIFYING");
        let now = Instant::now();
        system.verify_sumcheck_transposed(&pols, f(0)).unwrap();
        println!("Verification took {:?}", now.elapsed());
    }
}
