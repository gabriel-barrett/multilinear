use anyhow::{ensure, Context, Result};
use rand_chacha::ChaCha8Rng;

use crate::{
    constraints::BQCS,
    fields::{random_base, BaseField},
    polynomials::{SquarePolynomial, SquarePolynomialEval},
};

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
        let mut acc = BaseField::from(0);
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
        let mut acc = BaseField::from(0);
        for i in 0..1 << (I - len) {
            let d = self.delta.evaluate_hypercube(row1, row2, i);
            if d == BaseField::from(0) {
                continue;
            }
            for j1 in 0..1 << J {
                let w1 = self.trace.evaluate_hypercube_col(row1, i, j1);
                for j2 in 0..1 << J {
                    let w2 = self.trace.evaluate_hypercube_col(row2, i, j2);
                    let a = self.matrix.evaluate_hypercube(&[], j1, &[], j2);
                    acc += d * a * w1 * w2
                }
            }
        }
        acc
    }

    pub fn partial_sum_cols(&self, col1: &[BaseField], col2: &[BaseField]) -> BaseField {
        let len = col1.len();
        assert_eq!(col2.len(), len);
        assert!(len > 0);
        assert!(len <= J);
        let mut acc = BaseField::from(0);
        for i in 0..1 << I {
            let d = self.delta.evaluate_hypercube(&[], &[], i);
            if d == BaseField::from(0) {
                continue;
            }
            for j1 in 0..1 << (J - len) {
                let w1 = self.trace.evaluate_hypercube_row(i, col1, j1);
                for j2 in 0..1 << (J - len) {
                    let w2 = self.trace.evaluate_hypercube_row(i, col2, j2);
                    let a = self.matrix.evaluate_hypercube(col1, j1, col2, j2);
                    acc += d * a * w1 * w2;
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
        let mut acc = BaseField::from(0);
        let d = self.delta.evaluate_hypercube(row1, row2, 0);
        if d == BaseField::from(0) {
            return BaseField::from(0);
        }
        for j1 in 0..1 << (J - len) {
            let w1 = self.trace.evaluate_hypercube(row1, 0, col1, j1);
            for j2 in 0..1 << (J - len) {
                let a = self.matrix.evaluate_hypercube(col1, j1, col2, j2);
                let w2 = self.trace.evaluate_hypercube(row2, 0, col2, j2);
                acc += d * a * w1 * w2
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
        let mut acc = BaseField::from(0);
        let a = self.matrix.evaluate_hypercube(col1, 0, col2, 0);
        for i in 0..1 << (I - len) {
            let d = self.delta.evaluate_hypercube(row1, row2, i);
            if d == BaseField::from(0) {
                continue;
            }
            let w1 = self.trace.evaluate_hypercube(row1, i, col1, 0);
            let w2 = self.trace.evaluate_hypercube(row2, i, col2, 0);
            acc += d * a * w1 * w2
        }
        acc
    }

    pub fn generate_partial_polynomials(&mut self) -> Vec<SquarePolynomial<BaseField>> {
        let mut acc = Vec::with_capacity(2 * I + 2 * J);
        let mut r_acc = Vec::with_capacity(2 * I + 2 * J);

        let mut index_j1 = Vec::with_capacity(J);
        let mut index_j2 = Vec::with_capacity(J);
        for _ in 0..J {
            self.push_polynomials(
                &mut index_j1,
                &mut index_j2,
                &mut acc,
                &mut r_acc,
                |cs, index1, index2| cs.partial_sum_cols(index1, index2),
            );
        }

        let mut iter = r_acc.iter();
        let cols = [(); J].map(|_| (*iter.next().unwrap(), *iter.next().unwrap()));
        let col1 = &cols.map(|(col, _)| col);
        let col2 = &cols.map(|(_, col)| col);
        let mut index_i1 = Vec::with_capacity(I);
        let mut index_i2 = Vec::with_capacity(I);
        for _ in 0..I {
            self.push_polynomials(
                &mut index_i1,
                &mut index_i2,
                &mut acc,
                &mut r_acc,
                |cs, index1, index2| cs.partial_sum_rows_fixed_col(index1, index2, col1, col2),
            );
        }

        acc
    }

    pub fn generate_partial_polynomials_transposed(&mut self) -> Vec<SquarePolynomial<BaseField>> {
        let mut acc = Vec::with_capacity(2 * I + 2 * J);
        let mut r_acc = Vec::with_capacity(2 * I + 2 * J);

        let mut index_i1 = Vec::with_capacity(I);
        let mut index_i2 = Vec::with_capacity(I);
        for _ in 0..I {
            self.push_polynomials(
                &mut index_i1,
                &mut index_i2,
                &mut acc,
                &mut r_acc,
                |cs, index1, index2| cs.partial_sum_rows(index1, index2),
            );
        }

        let mut iter = r_acc.iter();
        let rows = [(); I].map(|_| (*iter.next().unwrap(), *iter.next().unwrap()));
        let row1 = &rows.map(|(row, _)| row);
        let row2 = &rows.map(|(_, row)| row);
        let mut index_j1 = Vec::with_capacity(J);
        let mut index_j2 = Vec::with_capacity(J);
        for _ in 0..J {
            self.push_polynomials(
                &mut index_j1,
                &mut index_j2,
                &mut acc,
                &mut r_acc,
                |cs, index1, index2| cs.partial_sum_cols_fixed_row(row1, row2, index1, index2),
            );
        }

        acc
    }

    fn push_polynomials(
        &mut self,
        index1: &mut Vec<BaseField>,
        index2: &mut Vec<BaseField>,
        acc: &mut Vec<SquarePolynomial<BaseField>>,
        r_acc: &mut Vec<BaseField>,
        partial_sum: impl Fn(&Self, &[BaseField], &[BaseField]) -> BaseField,
    ) {
        let f = BaseField::from;
        index1.push(f(-1));
        index2.push(f(0));
        let minus_one0 = partial_sum(self, index1, index2);
        let minus_one1 = partial_sum(self, index1, modify_last(index2, f(1)));
        let zero1 = partial_sum(self, modify_last(index1, f(0)), index2);
        let zero0 = partial_sum(self, index1, modify_last(index2, f(0)));
        let one0 = partial_sum(self, modify_last(index1, f(1)), index2);
        let one1 = partial_sum(self, index1, modify_last(index2, f(1)));
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
        let r1 = random_base(&mut self.rng);
        let pol2 = SquarePolynomialEval {
            minus_one: partial_sum(self, modify_last(index1, r1), modify_last(index2, f(-1))),
            zero: pol10.evaluate(r1),
            one: pol11.evaluate(r1),
        };
        let r2 = random_base(&mut self.rng);
        modify_last(index2, r2);
        acc.push(pol1.to_polynomial());
        acc.push(pol2.to_polynomial());
        r_acc.push(r1);
        r_acc.push(r2);
    }

    #[allow(dead_code)]
    fn push_polynomials_no_opt(
        &mut self,
        index1: &mut Vec<BaseField>,
        index2: &mut Vec<BaseField>,
        acc: &mut Vec<SquarePolynomial<BaseField>>,
        r_acc: &mut Vec<BaseField>,
        partial_sum: impl Fn(&Self, &[BaseField], &[BaseField]) -> BaseField,
    ) {
        let f = BaseField::from;
        index1.push(f(-1));
        index2.push(f(0));
        let minus_one = partial_sum(self, index1, index2)
            + partial_sum(self, index1, modify_last(index2, f(1)));
        let zero = partial_sum(self, modify_last(index1, f(0)), index2)
            + partial_sum(self, index1, modify_last(index2, f(0)));
        let one = partial_sum(self, modify_last(index1, f(1)), index2)
            + partial_sum(self, index1, modify_last(index2, f(1)));
        let pol1 = SquarePolynomialEval {
            minus_one,
            zero,
            one,
        };
        let r1 = random_base(&mut self.rng);
        let pol2 = SquarePolynomialEval {
            minus_one: partial_sum(self, modify_last(index1, r1), modify_last(index2, f(-1))),
            zero: partial_sum(self, index1, modify_last(index2, f(0))),
            one: partial_sum(self, index1, modify_last(index2, f(1))),
        };
        let r2 = random_base(&mut self.rng);
        modify_last(index2, r2);
        acc.push(pol1.to_polynomial());
        acc.push(pol2.to_polynomial());
        r_acc.push(r1);
        r_acc.push(r2);
    }

    pub fn verify_sumcheck(
        &self,
        rng: &mut ChaCha8Rng,
        pols: &[SquarePolynomial<BaseField>],
        sum: BaseField,
    ) -> Result<()> {
        let mut row1 = [BaseField::from(0); I];
        let mut row2 = [BaseField::from(0); I];
        let mut col1 = [BaseField::from(0); J];
        let mut col2 = [BaseField::from(0); J];
        let mut rs = Vec::with_capacity(2 * I + 2 * J);
        (0..J).for_each(|j| {
            let r1 = random_base(rng);
            let r2 = random_base(rng);
            col1[j] = r1;
            col2[j] = r2;
            rs.push(r1);
            rs.push(r2);
        });
        (0..I).for_each(|i| {
            let r1 = random_base(rng);
            let r2 = random_base(rng);
            row1[i] = r1;
            row2[i] = r2;
            rs.push(r1);
            rs.push(r2);
        });
        self.verify_sumcheck_constraints(pols, &rs, sum, &row1, &row2, &col1, &col2)
    }

    pub fn verify_sumcheck_transposed(
        &self,
        rng: &mut ChaCha8Rng,
        pols: &[SquarePolynomial<BaseField>],
        sum: BaseField,
    ) -> Result<()> {
        let mut row1 = [BaseField::from(0); I];
        let mut row2 = [BaseField::from(0); I];
        let mut col1 = [BaseField::from(0); J];
        let mut col2 = [BaseField::from(0); J];
        let mut rs = Vec::with_capacity(2 * I + 2 * J);
        (0..I).for_each(|i| {
            let r1 = random_base(rng);
            let r2 = random_base(rng);
            row1[i] = r1;
            row2[i] = r2;
            rs.push(r1);
            rs.push(r2);
        });
        (0..J).for_each(|j| {
            let r1 = random_base(rng);
            let r2 = random_base(rng);
            col1[j] = r1;
            col2[j] = r2;
            rs.push(r1);
            rs.push(r2);
        });
        self.verify_sumcheck_constraints(pols, &rs, sum, &row1, &row2, &col1, &col2)
    }

    #[allow(clippy::too_many_arguments)]
    fn verify_sumcheck_constraints(
        &self,
        pols: &[SquarePolynomial<BaseField>],
        rs: &[BaseField],
        sum: BaseField,
        row1: &[BaseField; I],
        row2: &[BaseField; I],
        col1: &[BaseField; J],
        col2: &[BaseField; J],
    ) -> Result<()> {
        assert_eq!(pols.len(), 2 * I + 2 * J);
        let f = BaseField::from;
        let mut iter = rs.iter().zip(pols.iter());
        let (mut r, mut pol) = iter.next().context("Expects at least one polynomial")?;
        ensure!(
            sum == pol.evaluate(f(0)) + pol.evaluate(f(1)),
            "Does not sumcheck"
        );
        for (i, (r_next, pol_next)) in iter.enumerate() {
            ensure!(
                pol.evaluate(*r) == pol_next.evaluate(f(0)) + pol_next.evaluate(f(1)),
                "Polynomial {i} does not satisfy equation"
            );
            r = r_next;
            pol = pol_next;
        }
        ensure!(
            self.evaluate(row1, row2, col1, col2) == pol.evaluate(*r),
            "Does not match polynomial evaluation"
        );
        Ok(())
    }
}

// Auxiliary function
fn modify_last<T>(vec: &mut [T], a: T) -> &mut [T] {
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
        let system = pythagorean_cs();
        assert_eq!(system.full_sum(), BaseField::from(0));
    }

    #[test]
    fn sumcheck_test() {
        let mut system = pythagorean_cs();
        let mut rng = system.rng.clone();
        let pols = system.generate_partial_polynomials();
        println!("{pols:#?}");
        system
            .verify_sumcheck(&mut rng, &pols, BaseField::from(0))
            .unwrap();
    }

    #[test]
    fn sumcheck_high_bench() {
        let mut rows = pythagorean_trace().rows.to_vec();
        const TOTAL_LOG_HEIGHT: usize = 18;
        for _ in 0..TOTAL_LOG_HEIGHT - I {
            rows.extend(rows.clone());
        }
        assert_eq!(rows.len(), 1 << TOTAL_LOG_HEIGHT);
        let trace = Trace::<TOTAL_LOG_HEIGHT, J> {
            rows: rows.try_into().unwrap(),
        };
        let set = pythagorean_set();
        let mut system = BQCS::new(trace, set);
        let mut rng = system.rng.clone();
        println!(
            "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
            1 << TOTAL_LOG_HEIGHT,
            1 << J,
        );
        let now = Instant::now();
        let pols = system.generate_partial_polynomials();
        println!("Generation took {:?}", now.elapsed());
        println!("VERIFYING");
        let now = Instant::now();
        system
            .verify_sumcheck(&mut rng, &pols, BaseField::from(0))
            .unwrap();
        println!("Verification took {:?}", now.elapsed());
    }

    #[test]
    fn sumcheck_high_transposed_bench() {
        let mut rows = pythagorean_trace().rows.to_vec();
        const TOTAL_LOG_HEIGHT: usize = 18;
        for _ in 0..TOTAL_LOG_HEIGHT - I {
            rows.extend(rows.clone());
        }
        assert_eq!(rows.len(), 1 << TOTAL_LOG_HEIGHT);
        let trace = Trace::<TOTAL_LOG_HEIGHT, J> {
            rows: rows.try_into().unwrap(),
        };
        let set = pythagorean_set();
        let mut system = BQCS::new(trace, set);
        let mut rng = system.rng.clone();
        println!(
            "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
            1 << TOTAL_LOG_HEIGHT,
            1 << J,
        );
        let now = Instant::now();
        let pols = system.generate_partial_polynomials_transposed();
        println!("Generation took {:?}", now.elapsed());
        println!("VERIFYING");
        let now = Instant::now();
        system
            .verify_sumcheck_transposed(&mut rng, &pols, BaseField::from(0))
            .unwrap();
        println!("Verification took {:?}", now.elapsed());
    }
}
