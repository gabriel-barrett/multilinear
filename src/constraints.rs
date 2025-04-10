use crate::BaseField;

#[derive(Clone, Copy)]
pub struct Var<const J: usize>(usize);

impl<const J: usize> Var<J> {
    pub fn new(value: usize) -> Self {
        assert!(value < 1 << J);
        Self(value)
    }

    pub fn evaluate(&self, points: &[BaseField; J]) -> BaseField {
        let index = self.0;
        let one = BaseField::from(1);
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

pub struct LinearCombination<const J: usize>(Box<[(BaseField, Var<J>)]>);

impl<const J: usize> LinearCombination<J> {
    pub fn evaluate(&self, points: &[BaseField; J]) -> BaseField {
        self.0.iter().map(|(c, x)| *c * x.evaluate(points)).sum()
    }
}

pub struct Constraint<const J: usize>(Box<[(LinearCombination<J>, LinearCombination<J>)]>);

impl<const J: usize> Constraint<J> {
    pub fn evaluate(&self, col1: &[BaseField; J], col2: &[BaseField; J]) -> BaseField {
        self.0
            .iter()
            .map(|(lc1, lc2)| lc1.evaluate(col1) * lc2.evaluate(col2))
            .sum()
    }
}

pub struct ConstraintSet<const J: usize, const K: usize>
where
    [(); 1 << K]:,
{
    constraints: Box<[Constraint<J>; 1 << K]>,
}

impl<const J: usize, const K: usize> ConstraintSet<J, K>
where
    [(); 1 << K]:,
{
    pub fn evaluate(
        &self,
        col1: &[BaseField; J],
        col2: &[BaseField; J],
        constraint: &[BaseField; K],
    ) -> BaseField {
        self.constraints
            .iter()
            .enumerate()
            .map(|(k, expression)| {
                let constraint_mask = Var::<K>(k).evaluate(constraint);
                constraint_mask * expression.evaluate(col1, col2)
            })
            .sum()
    }
}

pub struct Trace<const I: usize, const J: usize>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    rows: Box<[[BaseField; 1 << J]; 1 << I]>,
}

impl<const I: usize, const J: usize> Trace<I, J>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    pub fn evaluate(&self, row: &[BaseField; I], col: &[BaseField; J]) -> BaseField {
        let mut res = BaseField::from(0);
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

pub struct Delta<const I: usize> {
    data: [BaseField; I],
}

impl<const I: usize> Delta<I> {
    pub fn evaluate_unary(&self, b: &[BaseField; I]) -> BaseField {
        let one = BaseField::from(1);
        let pass = |i| {
            let a = self.data[i];
            let b = b[i];
            a * b + (one - a) * (one - b)
        };
        (0..I).map(pass).product()
    }

    pub fn evaluate_binary(&self, b: &[BaseField; I], c: &[BaseField; I]) -> BaseField {
        let one = BaseField::from(1);
        let pass = |i| {
            let a = self.data[i];
            let b = b[i];
            let c = c[i];
            a * b * c + (one - a) * (one - b) * (one - c)
        };
        (0..I).map(pass).product()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trace_evaluation_test() {
        let f = BaseField::from;
        let rows = [
            [1, 2, 3, 4].map(f),
            [5, 6, 7, 8].map(f),
            [9, 10, 11, 12].map(f),
            [13, 14, 15, 16].map(f),
        ]
        .into();
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
        assert_eq!(delta.evaluate_binary(&point1, &point2), f(0));
        let res = f((1 - 4) * 8 * 12 * (1 - 10));
        assert_eq!(delta.evaluate_binary(&point1, &point1), res);
        assert_eq!(delta.evaluate_unary(&point1), res);
    }
}
