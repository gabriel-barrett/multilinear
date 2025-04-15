#![allow(non_snake_case)]
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
    pub fn generate_partial_polynomials(
        &mut self,
        sum: BaseField,
    ) -> Vec<SquarePolynomial<BaseField>> {
        let mut delta = self.delta.build_table();
        let mut matrix = self.matrix.build_table();
        let mut trace1 = self.trace.build_table();
        let mut trace2 = self.trace.build_table();
        let mut pols = Vec::with_capacity(2 * I + 2 * J);
        let mut previous_sum = sum;
        for _ in 0..J {
            let a00 = partial_sum_cols(0, 0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a01 = partial_sum_cols(0, 1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a10 = partial_sum_cols(1, 0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a11 = previous_sum - a00 - a01 - a10;
            let b0 = compute_minus_one_cols(0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let b1 = compute_minus_one_cols(1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let pol1 = SquarePolynomialEval {
                zero: a00 + a01,
                one: a10 + a11,
                minus_one: b0 + b1,
            };
            let r1 = self.random();
            let minus_one =
                compute_r_minus_one_cols(r1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let zero = SquarePolynomialEval {
                zero: a00,
                one: a10,
                minus_one: b0,
            }
            .evaluate(r1);
            let one = SquarePolynomialEval {
                zero: a01,
                one: a11,
                minus_one: b1,
            }
            .evaluate(r1);
            let pol2 = SquarePolynomialEval {
                zero,
                one,
                minus_one,
            };
            let r2 = self.random();
            pols.push(pol1.to_polynomial());
            pols.push(pol2.to_polynomial());
            previous_sum = pol2.evaluate(r2);
            fold_matrices_col(r1, r2, &mut delta, &mut matrix, &mut trace1, &mut trace2);
        }
        for _ in 0..I {
            let a00 = partial_sum_rows(0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a01 = BaseField::from(0);
            let a10 = BaseField::from(0);
            let a11 = previous_sum - a00 - a01 - a10;
            let b0 = compute_minus_one_rows(0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let b1 = compute_minus_one_rows(1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let pol1 = SquarePolynomialEval {
                zero: a00 + a01,
                one: a10 + a11,
                minus_one: b0 + b1,
            };
            let r1 = self.random();
            let minus_one =
                compute_r_minus_one_rows(r1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let zero = SquarePolynomialEval {
                zero: a00,
                one: a10,
                minus_one: b0,
            }
            .evaluate(r1);
            let one = SquarePolynomialEval {
                zero: a01,
                one: a11,
                minus_one: b1,
            }
            .evaluate(r1);
            let pol2 = SquarePolynomialEval {
                zero,
                one,
                minus_one,
            };
            let r2 = self.random();
            pols.push(pol1.to_polynomial());
            pols.push(pol2.to_polynomial());
            previous_sum = pol2.evaluate(r2);
            fold_matrices_row(r1, r2, &mut delta, &mut matrix, &mut trace1, &mut trace2);
        }
        pols
    }

    pub fn generate_partial_polynomials_transposed(
        &mut self,
        sum: BaseField,
    ) -> Vec<SquarePolynomial<BaseField>> {
        let mut delta = self.delta.build_table();
        let mut matrix = self.matrix.build_table();
        let mut trace1 = self.trace.build_table();
        let mut trace2 = self.trace.build_table();
        let mut pols = Vec::with_capacity(2 * I + 2 * J);
        let mut previous_sum = sum;
        for _ in 0..I {
            let a00 = partial_sum_rows(0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a01 = BaseField::from(0);
            let a10 = BaseField::from(0);
            let a11 = previous_sum - a00 - a01 - a10;
            let b0 = compute_minus_one_rows(0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let b1 = compute_minus_one_rows(1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let pol1 = SquarePolynomialEval {
                zero: a00 + a01,
                one: a10 + a11,
                minus_one: b0 + b1,
            };
            let r1 = self.random();
            let minus_one =
                compute_r_minus_one_rows(r1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let zero = SquarePolynomialEval {
                zero: a00,
                one: a10,
                minus_one: b0,
            }
            .evaluate(r1);
            let one = SquarePolynomialEval {
                zero: a01,
                one: a11,
                minus_one: b1,
            }
            .evaluate(r1);
            let pol2 = SquarePolynomialEval {
                zero,
                one,
                minus_one,
            };
            let r2 = self.random();
            pols.push(pol1.to_polynomial());
            pols.push(pol2.to_polynomial());
            previous_sum = pol2.evaluate(r2);
            fold_matrices_row(r1, r2, &mut delta, &mut matrix, &mut trace1, &mut trace2);
        }
        for _ in 0..J {
            let a00 = partial_sum_cols(0, 0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a01 = partial_sum_cols(0, 1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a10 = partial_sum_cols(1, 0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let a11 = previous_sum - a00 - a01 - a10;
            let b0 = compute_minus_one_cols(0, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let b1 = compute_minus_one_cols(1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let pol1 = SquarePolynomialEval {
                zero: a00 + a01,
                one: a10 + a11,
                minus_one: b0 + b1,
            };
            let r1 = self.random();
            let minus_one =
                compute_r_minus_one_cols(r1, &mut delta, &mut matrix, &mut trace1, &mut trace2);
            let zero = SquarePolynomialEval {
                zero: a00,
                one: a10,
                minus_one: b0,
            }
            .evaluate(r1);
            let one = SquarePolynomialEval {
                zero: a01,
                one: a11,
                minus_one: b1,
            }
            .evaluate(r1);
            let pol2 = SquarePolynomialEval {
                zero,
                one,
                minus_one,
            };
            let r2 = self.random();
            pols.push(pol1.to_polynomial());
            pols.push(pol2.to_polynomial());
            previous_sum = pol2.evaluate(r2);
            fold_matrices_col(r1, r2, &mut delta, &mut matrix, &mut trace1, &mut trace2);
        }
        pols
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

fn compute_minus_one_rows(
    i_bit: usize,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) -> BaseField {
    let I = delta.len() >> 1;
    let J = matrix.len();
    let i_offset = 1 << I.trailing_zeros();
    let maybe_offset = i_offset * i_bit;

    let mut acc = BaseField::from(0);
    let two = BaseField::from(2);
    let minus_one = BaseField::from(-1);
    let delta_coeff = if maybe_offset == 0 {
        two
    } else {
        minus_one
    };
    for i in 0..I {
        let d = delta[i + maybe_offset];
        for j1 in 0..J {
            let w1 = two * trace1[i][j1] - trace1[i + i_offset][j1];
            for j2 in 0..J {
                let w2 = trace2[i + maybe_offset][j2];
                let a = matrix[j1][j2];
                acc += d * a * w1 * w2;
            }
        }
    }
    delta_coeff * acc
}

fn compute_r_minus_one_rows(
    r1: BaseField,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) -> BaseField {
    let I = delta.len() >> 1;
    let J = matrix.len();
    let i_offset = 1 << I.trailing_zeros();

    let two = BaseField::from(2);
    let s1 = BaseField::from(1) - r1;
    let coeff1 = two * s1;
    let coeff2 = r1;
    let mut acc = BaseField::from(0);
    for i in 0..I {
        let d = coeff1 * delta[i] - coeff2 * delta[i + i_offset];
        for j1 in 0..J {
            let w1 = s1 * trace1[i][j1] + r1 * trace1[i + i_offset][j1];
            for j2 in 0..0 + J {
                let w2 = two * trace2[i][j2] - trace2[i + i_offset][j2];
                let a = matrix[j1][j2];
                acc += d * a * w1 * w2;
            }
        }
    }
    acc
}

fn compute_minus_one_cols(
    j2_bit: usize,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) -> BaseField {
    let I = delta.len();
    let J = matrix.len() >> 1;
    let j1_offset = 1 << J.trailing_zeros();
    let j2_offset = j2_bit << J.trailing_zeros();

    let two = BaseField::from(2);
    let mut acc = BaseField::from(0);
    for i in 0..I {
        let d = delta[i];
        for j1 in 0..J {
            let w1_0 = trace1[i][j1];
            let w1_1 = trace1[i][j1 + j1_offset];
            let w1 = two * w1_0 - w1_1;
            for j2 in j2_offset..j2_offset + J {
                let w2 = trace2[i][j2];
                let a_0 = matrix[j1][j2];
                let a_1 = matrix[j1 + j1_offset][j2];
                let a = two * a_0 - a_1;
                acc += d * a * w1 * w2;
            }
        }
    }
    acc
}

fn compute_r_minus_one_cols(
    r1: BaseField,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) -> BaseField {
    let I = delta.len();
    let J = matrix.len() >> 1;
    let offset = 1 << J.trailing_zeros();

    let two = BaseField::from(2);
    let s1 = BaseField::from(1) - r1;
    let mut acc = BaseField::from(0);
    let coeff1 = s1 * two;
    let coeff2 = s1;
    let coeff3 = r1 * two;
    let coeff4 = r1;
    for i in 0..I {
        let d = delta[i];
        for j1 in 0..J {
            let w1_0 = s1 * trace1[i][j1];
            let w1_1 = r1 * trace1[i][j1 + offset];
            let w1 = w1_0 + w1_1;
            for j2 in 0..J {
                let w2_0 = two * trace2[i][j2];
                let w2_1 = trace2[i][j2 + offset];
                let w2 = w2_0 - w2_1;
                let a_1 = coeff1 * matrix[j1][j2];
                let a_2 = coeff2 * matrix[j1][j2 + offset];
                let a_3 = coeff3 * matrix[j1 + offset][j2];
                let a_4 = coeff4 * matrix[j1 + offset][j2 + offset];
                let a = a_1 - a_2 + a_3 - a_4;
                acc += d * a * w1 * w2;
            }
        }
    }
    acc
}

fn partial_sum_cols(
    j1_bit: usize,
    j2_bit: usize,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) -> BaseField {
    let I = delta.len();
    let J = matrix.len() >> 1;
    let j1_offset = j1_bit << J.trailing_zeros();
    let j2_offset = j2_bit << J.trailing_zeros();

    let mut acc = BaseField::from(0);
    for i in 0..I {
        let d = delta[i];
        for j1 in j1_offset..j1_offset + J {
            let w1 = trace1[i][j1];
            for j2 in j2_offset..j2_offset + J {
                let w2 = trace2[i][j2];
                let a = matrix[j1][j2];
                acc += d * a * w1 * w2;
            }
        }
    }
    acc
}

fn partial_sum_rows(
    i_bit: usize,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) -> BaseField {
    let I = delta.len() >> 1;
    let J = matrix.len();
    let i_offset = i_bit << I.trailing_zeros();

    let mut acc = BaseField::from(0);
    for i in i_offset..i_offset + I {
        let d = delta[i];
        for j1 in 0..J {
            let w1 = trace1[i][j1];
            for j2 in 0..J {
                let w2 = trace2[i][j2];
                let a = matrix[j1][j2];
                acc += d * a * w1 * w2;
            }
        }
    }
    acc
}

fn fold_matrices_row(
    r1: BaseField,
    r2: BaseField,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) {
    let I = delta.len() >> 1;
    let J = matrix.len();
    let i_offset = 1 << I.trailing_zeros();
    let s1 = BaseField::from(1) - r1;
    let s2 = BaseField::from(1) - r2;
    let coeff1 = s1 * s2;
    let coeff4 = r1 * r2;
    for i in 0..I {
        delta[i] = coeff1 * delta[i] + coeff4 * delta[i + i_offset];
        for j1 in 0..J {
            trace1[i][j1] = s1 * trace1[i][j1] + r1 * trace1[i + i_offset][j1];
        }
        for j2 in 0..J {
            trace2[i][j2] = s2 * trace2[i][j2] + r2 * trace2[i + i_offset][j2];
        }
    }
    delta.truncate(I);
    trace1.truncate(I);
    trace2.truncate(I);
}

fn fold_matrices_col(
    r1: BaseField,
    r2: BaseField,
    delta: &mut Vec<BaseField>,
    matrix: &mut Vec<Vec<BaseField>>,
    trace1: &mut Vec<Vec<BaseField>>,
    trace2: &mut Vec<Vec<BaseField>>,
) {
    let I = delta.len();
    let J = matrix.len() >> 1;
    let j_offset = 1 << J.trailing_zeros();
    let s1 = BaseField::from(1) - r1;
    let s2 = BaseField::from(1) - r2;
    let coeff1 = s1 * s2;
    let coeff2 = r1 * s2;
    let coeff3 = s1 * r2;
    let coeff4 = r1 * r2;
    for i in 0..I {
        for j1 in 0..J {
            trace1[i][j1] = s1 * trace1[i][j1] + r1 * trace1[i][j1 + j_offset];
        }
        trace1[i].truncate(J);
        for j2 in 0..J {
            trace2[i][j2] = s2 * trace2[i][j2] + r2 * trace2[i][j2 + j_offset];
        }
        trace2[i].truncate(J);
    }
    for j1 in 0..J {
        for j2 in 0..J {
            matrix[j1][j2] = coeff1 * matrix[j1][j2]
                + coeff2 * matrix[j1 + j_offset][j2]
                + coeff3 * matrix[j1][j2 + j_offset]
                + coeff4 * matrix[j1 + j_offset][j2 + j_offset];
        }
        matrix[j1].truncate(J);
    }
    matrix.truncate(J);
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
    fn sumcheck_test() {
        let mut system = pythagorean_cs();
        let mut rng = system.rng.clone();
        let sum = BaseField::from(0);
        let pols = system.generate_partial_polynomials(sum);
        println!("{pols:#?}");
        system.verify_sumcheck(&mut rng, &pols, sum).unwrap();
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
        let sum = BaseField::from(0);
        println!(
            "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
            1 << TOTAL_LOG_HEIGHT,
            1 << J,
        );
        let now = Instant::now();
        let pols = system.generate_partial_polynomials(sum);
        println!("Generation took {:?}", now.elapsed());
        println!("VERIFYING");
        let now = Instant::now();
        system.verify_sumcheck(&mut rng, &pols, sum).unwrap();
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
        let sum = BaseField::from(0);
        println!(
            "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
            1 << TOTAL_LOG_HEIGHT,
            1 << J,
        );
        let now = Instant::now();
        let pols = system.generate_partial_polynomials_transposed(sum);
        println!("Generation took {:?}", now.elapsed());
        println!("VERIFYING");
        let now = Instant::now();
        system
            .verify_sumcheck_transposed(&mut rng, &pols, sum)
            .unwrap();
        println!("Verification took {:?}", now.elapsed());
    }
}
