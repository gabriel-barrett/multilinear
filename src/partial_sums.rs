use anyhow::{ensure, Context, Result};
use std::ops::Add;

use crate::{
    constraints::{ConstraintSet, ConstraintSetMatrix, Delta, Trace},
    BaseField,
};
use rand::Rng;

// Batched quadratic constraint system
pub struct BQCS<const I: usize, const J: usize, const K: usize>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
    [(); 1 << K]:,
{
    pub delta: Delta<I>,
    pub trace: Trace<I, J>,
    pub matrix: ConstraintSetMatrix<J, K>,
}

impl<const I: usize, const J: usize, const K: usize> BQCS<I, J, K>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
    [(); 1 << K]:,
{
    pub fn new(trace: Trace<I, J>, set: ConstraintSet<J, K>) -> Self {
        let (random_i, random_k) = generate_challenges(&trace);
        let delta = Delta { data: random_i };
        let matrix = ConstraintSetMatrix { set, random_k };
        Self {
            delta,
            trace,
            matrix,
        }
    }

    pub fn evaluate(
        &self,
        row1: &[BaseField; I],
        row2: &[BaseField; I],
        col1: &[BaseField; J],
        col2: &[BaseField; J],
    ) -> BaseField {
        let d = self.delta.evaluate_binary(row1, row2);
        let a = self.matrix.evaluate(col1, col2);
        let w1 = self.trace.evaluate(row1, col1);
        let w2 = self.trace.evaluate(row2, col2);
        d * a * w1 * w2
    }

    pub fn evaluate_unary(
        &self,
        row: &[BaseField; I],
        col1: &[BaseField; J],
        col2: &[BaseField; J],
    ) -> BaseField {
        let d = self.delta.evaluate_unary(row);
        let a = self.matrix.evaluate(col1, col2);
        let w1 = self.trace.evaluate(row, col1);
        let w2 = self.trace.evaluate(row, col2);
        d * a * w1 * w2
    }

    // Okay as it is linear in the height, but it is still width squared,
    // even if the constraint matrix usually has a sparse representation.
    // Needs further optimizations.
    pub fn full_sum(&self) -> BaseField {
        let mut acc = BaseField::from(0);
        for i in 0..1 << I {
            let i = to_hypercube_arr(i);
            for j1 in 0..1 << J {
                let j1 = to_hypercube_arr(j1);
                for j2 in 0..1 << J {
                    let j2 = to_hypercube_arr(j2);
                    acc += self.evaluate_unary(&i, &j1, &j2);
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
            let mut index_i1 = row1.to_vec();
            index_i1.extend(to_hypercube(i, I - len));
            let index_i1 = (&index_i1[..]).try_into().unwrap();
            let mut index_i2 = row2.to_vec();
            index_i2.extend(to_hypercube(i, I - len));
            let index_i2 = (&index_i2[..]).try_into().unwrap();
            for j1 in 0..1 << J {
                let j1 = to_hypercube_arr(j1);
                for j2 in 0..1 << J {
                    let j2 = to_hypercube_arr(j2);
                    acc += self.evaluate(index_i1, index_i2, &j1, &j2);
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
        for j1 in 0..1 << (J - len) {
            let mut index_j1 = col1.to_vec();
            index_j1.extend(to_hypercube(j1, J - len));
            let index_j1 = (&index_j1[..]).try_into().unwrap();
            for j2 in 0..1 << (J - len) {
                let mut index_j2 = col2.to_vec();
                index_j2.extend(to_hypercube(j2, J - len));
                let index_j2 = (&index_j2[..]).try_into().unwrap();
                acc += self.evaluate(row1, row2, &index_j1, &index_j2);
            }
        }
        acc
    }

    pub fn generate_partial_polynomials(&self) -> Vec<(BaseField, SquarePolynomial)> {
        fn modify_last(vec: &mut [BaseField], a: BaseField) {
            let len = vec.len();
            vec[len - 1] = a;
        }

        let f = BaseField::from;
        let mut rng = rand::thread_rng();
        let mut index_i1 = Vec::with_capacity(I);
        let mut index_i2 = Vec::with_capacity(I);
        let mut acc = Vec::with_capacity(2 * I + 2 * J);
        let mut row1 = Vec::with_capacity(I);
        let mut row2 = Vec::with_capacity(I);
        for _ in 0..I {
            index_i1.push(f(-1));
            index_i2.push(f(0));
            let minus_one0 = self.partial_sum_rows(&index_i1, &index_i2);
            modify_last(&mut index_i2, f(1));
            let minus_one1 = self.partial_sum_rows(&index_i1, &index_i2);
            modify_last(&mut index_i1, f(0));
            let zero1 = self.partial_sum_rows(&index_i1, &index_i2);
            modify_last(&mut index_i2, f(0));
            let zero0 = self.partial_sum_rows(&index_i1, &index_i2);
            modify_last(&mut index_i1, f(1));
            let one0 = self.partial_sum_rows(&index_i1, &index_i2);
            modify_last(&mut index_i2, f(1));
            let one1 = self.partial_sum_rows(&index_i1, &index_i2);
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
            let r1 = BaseField::from(rng.gen::<u64>());
            modify_last(&mut index_i1, r1);
            modify_last(&mut index_i2, f(-1));
            let pol2 = SquarePolynomialEval {
                minus_one: self.partial_sum_rows(&index_i1, &index_i2),
                zero: pol10.evaluate(r1),
                one: pol11.evaluate(r1),
            };
            let r2 = BaseField::from(rng.gen::<u64>());
            modify_last(&mut index_i2, r2);
            acc.push((r1, pol1.to_polynomial()));
            acc.push((r2, pol2.to_polynomial()));
            row1.push(r1);
            row2.push(r2);
        }
        let row1 = (&row1[..]).try_into().unwrap();
        let row2 = (&row2[..]).try_into().unwrap();
        let mut index_j1 = Vec::with_capacity(J);
        let mut index_j2 = Vec::with_capacity(J);
        for _ in 0..J {
            index_j1.push(f(-1));
            index_j2.push(f(0));
            let minus_one0 = self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2);
            modify_last(&mut index_j2, f(1));
            let minus_one1 = self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2);
            modify_last(&mut index_j1, f(0));
            let zero1 = self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2);
            modify_last(&mut index_j2, f(0));
            let zero0 = self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2);
            modify_last(&mut index_j1, f(1));
            let one0 = self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2);
            modify_last(&mut index_j2, f(1));
            let one1 = self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2);
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
            let r1 = BaseField::from(rng.gen::<u64>());
            modify_last(&mut index_j1, r1);
            modify_last(&mut index_j2, f(-1));
            let pol2 = SquarePolynomialEval {
                minus_one: self.partial_sum_cols_fixed_row(row1, row2, &index_j1, &index_j2),
                zero: pol10.evaluate(r1),
                one: pol11.evaluate(r1),
            };
            let r2 = BaseField::from(rng.gen::<u64>());
            modify_last(&mut index_j2, r2);
            acc.push((r1, pol1.to_polynomial()));
            acc.push((r2, pol2.to_polynomial()));
        }
        acc
    }

    pub fn verify_sumcheck(
        &self,
        pols: &[(BaseField, SquarePolynomial)],
        sum: BaseField,
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
        let mut r_iter = pols.iter().map(|(r, _)| r);
        let rows = [(); I].map(|_| (*r_iter.next().unwrap(), *r_iter.next().unwrap()));
        let cols = [(); J].map(|_| (*r_iter.next().unwrap(), *r_iter.next().unwrap()));
        let row1 = rows.map(|(row, _)| row);
        let row2 = rows.map(|(_, row)| row);
        let col1 = cols.map(|(col, _)| col);
        let col2 = cols.map(|(_, col)| col);
        ensure!(
            self.evaluate(&row1, &row2, &col1, &col2) == pol.evaluate(r),
            "Does not match polynomial evaluation"
        );
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SquarePolynomialEval {
    pub one: BaseField,
    pub zero: BaseField,
    pub minus_one: BaseField,
}

impl Add for SquarePolynomialEval {
    type Output = SquarePolynomialEval;
    fn add(self, other: SquarePolynomialEval) -> SquarePolynomialEval {
        SquarePolynomialEval {
            one: self.one + other.one,
            zero: self.zero + other.zero,
            minus_one: self.minus_one + other.minus_one,
        }
    }
}

impl SquarePolynomialEval {
    pub fn to_polynomial(self) -> SquarePolynomial {
        let a = self.zero;
        let b = (self.one - self.minus_one) / BaseField::from(2);
        let c = (self.one + self.minus_one) / BaseField::from(2) - a;
        SquarePolynomial { coeffs: [a, b, c] }
    }

    pub fn evaluate(&self, x: BaseField) -> BaseField {
        self.to_polynomial().evaluate(x)
    }
}

#[derive(Clone, Copy)]
pub struct SquarePolynomial {
    pub coeffs: [BaseField; 3],
}

impl std::fmt::Debug for SquarePolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [a, b, c] = self.coeffs;
        write!(f, "{} + {}*X + {}*X^2", a, b, c)?;
        Ok(())
    }
}

impl SquarePolynomial {
    pub fn evaluate(&self, x: BaseField) -> BaseField {
        let [a, b, c] = self.coeffs;
        a + b * x + c * x * x
    }
}

fn to_hypercube(index: u64, len: usize) -> Vec<BaseField> {
    let mut points = vec![BaseField::from(0); len];
    (0..len).for_each(|i| {
        if (index >> i) & 1 == 1 {
            points[len - 1 - i] = BaseField::from(1);
        }
    });
    points
}

fn to_hypercube_arr<const I: usize>(index: u64) -> [BaseField; I] {
    let mut points = [BaseField::from(0); I];
    (0..I).for_each(|i| {
        if (index >> i) & 1 == 1 {
            points[I - 1 - i] = BaseField::from(1);
        }
    });
    points
}

fn generate_challenges<const I: usize, const J: usize, const K: usize>(
    _trace: &Trace<I, J>,
) -> ([BaseField; I], [BaseField; K])
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    let mut rng = rand::thread_rng();
    let random_i = [(); I].map(|_| BaseField::from(rng.gen::<u64>()));
    let random_k = [(); K].map(|_| BaseField::from(rng.gen::<u64>()));
    (random_i, random_k)
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use crate::expr::Expr;

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
}
