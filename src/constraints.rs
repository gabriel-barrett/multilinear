use crate::BaseField;
use ark_ff::{AdditiveGroup, Field};
use rand::Rng;

#[derive(Clone, Copy, Debug)]
pub struct Var<const J: usize>(pub(crate) usize);

impl<const J: usize> Var<J> {
    pub fn new(value: usize) -> Self {
        assert!(value < 1 << J);
        Self(value)
    }

    pub fn evaluate(&self, points: &[BaseField; J]) -> BaseField {
        let index = self.0;
        let one = BaseField::ONE;
        let select = |i| {
            // Note: the points are read from last to first, since WHIR
            // is big endian and we want to follow the same convention
            let point = points[J - 1 - i];
            if (index >> i) & 1 == 1 {
                point
            } else {
                one - point
            }
        };
        (0..J).map(select).product()
    }
}

#[derive(Clone, Debug)]
pub struct LinearCombination<const J: usize>(pub(crate) Box<[(BaseField, Var<J>)]>);

impl<const J: usize> LinearCombination<J> {
    pub fn evaluate(&self, points: &[BaseField; J]) -> BaseField {
        self.0.iter().map(|(c, x)| *c * x.evaluate(points)).sum()
    }
}

#[derive(Clone, Debug)]
pub struct Quadratic<const J: usize>(
    pub(crate) Box<[(LinearCombination<J>, LinearCombination<J>)]>,
);

impl<const J: usize> Quadratic<J> {
    pub fn evaluate(&self, col1: &[BaseField; J], col2: &[BaseField; J]) -> BaseField {
        self.0
            .iter()
            .map(|(lc1, lc2)| lc1.evaluate(col1) * lc2.evaluate(col2))
            .sum()
    }
}

#[derive(Clone, Debug)]
pub struct ConstraintSet<const J: usize, const K: usize>
where
    [(); 1 << K]:,
{
    pub expressions: Box<[Quadratic<J>; 1 << K]>,
}

#[derive(Clone, Debug)]
pub struct ConstraintSetMatrix<const J: usize, const K: usize>
where
    [(); 1 << K]:,
{
    pub set: ConstraintSet<J, K>,
    pub random_k: [BaseField; K],
}

impl<const J: usize, const K: usize> ConstraintSetMatrix<J, K>
where
    [(); 1 << K]:,
{
    pub fn evaluate(&self, col1: &[BaseField; J], col2: &[BaseField; J]) -> BaseField {
        let random_k = &self.random_k;
        self.set
            .expressions
            .iter()
            .enumerate()
            .map(|(k, expression)| {
                let constraint_mask = Var::<K>(k).evaluate(random_k);
                constraint_mask * expression.evaluate(col1, col2)
            })
            .sum()
    }
}

#[derive(Clone, Debug)]
pub struct Trace<const I: usize, const J: usize>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    pub rows: Box<[[BaseField; 1 << J]; 1 << I]>,
}

impl<const I: usize, const J: usize> Trace<I, J>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    pub fn evaluate(&self, row: &[BaseField; I], col: &[BaseField; J]) -> BaseField {
        let mut res = BaseField::ZERO;
        for (i, coeffs) in self.rows.iter().enumerate() {
            let row_mask = Var::<I>(i).evaluate(row);
            for (j, coeff) in coeffs.iter().enumerate() {
                let col_mask = Var::<J>(j).evaluate(col);
                res += *coeff * row_mask * col_mask;
            }
        }
        res
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Delta<const I: usize> {
    pub data: [BaseField; I],
}

impl<const I: usize> Delta<I> {
    pub fn evaluate(&self, b: &[BaseField; I], c: &[BaseField; I]) -> BaseField {
        let one = BaseField::ONE;
        let pass = |i| {
            let a = self.data[i];
            let b = b[i];
            let c = c[i];
            a * b * c + (one - a) * (one - b) * (one - c)
        };
        (0..I).map(pass).product()
    }
}

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
        let d = self.delta.evaluate(row1, row2);
        let a = self.matrix.evaluate(col1, col2);
        let c = d * a;
        if c == BaseField::ZERO {
            return BaseField::ZERO;
        }
        let w1 = self.trace.evaluate(row1, col1);
        let w2 = self.trace.evaluate(row2, col2);
        c * w1 * w2
    }
}

fn generate_challenges<const I: usize, const J: usize, const K: usize>(
    _trace: &Trace<I, J>,
) -> ([BaseField; I], [BaseField; K])
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    let mut rng = rand::rng();
    let random_i = [(); I].map(|_| BaseField::from(rng.random::<u64>()));
    let random_k = [(); K].map(|_| BaseField::from(rng.random::<u64>()));
    (random_i, random_k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trace_evaluation_test() {
        let f = BaseField::from;
        let rows = Box::new([
            [1, 2, 3, 4].map(f),
            [5, 6, 7, 8].map(f),
            [9, 10, 11, 12].map(f),
            [13, 14, 15, 16].map(f),
        ]);
        let trace = Trace { rows };
        let zero = [0, 0].map(f);
        let one = [0, 1].map(f);
        let two = [1, 0].map(f);
        let three = [1, 1].map(f);
        assert_eq!(trace.evaluate(&zero, &two), f(3));
        assert_eq!(trace.evaluate(&one, &one), f(6));
        assert_eq!(trace.evaluate(&two, &zero), f(9));
        assert_eq!(trace.evaluate(&three, &three), f(16));
    }

    #[test]
    fn delta_evaluation_test() {
        let f = BaseField::from;
        let delta = Delta {
            data: [4, 8, 12, 10].map(f),
        };
        let point1 = [0, 1, 1, 0].map(f);
        let point2 = [1, 1, 1, 0].map(f);
        assert_eq!(delta.evaluate(&point1, &point2), f(0));
        let res = f((1 - 4) * 8 * 12 * (1 - 10));
        assert_eq!(delta.evaluate(&point1, &point1), res);
    }
}
