use crate::BaseField;

#[derive(Clone, Copy, Debug)]
pub struct Var<const J: usize>(pub(crate) usize);

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

    // Assumes `bits` has J bits most, and that `J` < `usize::BITS`
    pub fn evaluate_hypercube(&self, points: &[BaseField], bits: usize) -> BaseField {
        let mask = (1 << (J - points.len())) - 1;
        let index = self.0;
        if bits != index & mask {
            return BaseField::from(0);
        }
        let one = BaseField::from(1);
        let mut acc = one;
        for (i, point) in points.iter().enumerate() {
            // Note: the points are read from last to first, since WHIR
            // is big endian and we want to follow the same convention
            if (index >> (J - 1 - i)) & 1 == 1 {
                acc *= point;
            } else {
                acc *= one - point;
            }
        }
        acc
    }
}

#[derive(Clone, Debug)]
pub struct LinearCombination<const J: usize>(pub(crate) Box<[(BaseField, Var<J>)]>);

impl<const J: usize> LinearCombination<J> {
    pub fn evaluate(&self, points: &[BaseField; J]) -> BaseField {
        self.0.iter().map(|(c, x)| *c * x.evaluate(points)).sum()
    }

    pub fn evaluate_hypercube(&self, points: &[BaseField], bits: usize) -> BaseField {
        self.0
            .iter()
            .map(|(c, x)| *c * x.evaluate_hypercube(points, bits))
            .sum()
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

    pub fn evaluate_hypercube(
        &self,
        col1: &[BaseField],
        bits1: usize,
        col2: &[BaseField],
        bits2: usize,
    ) -> BaseField {
        self.0
            .iter()
            .map(|(lc1, lc2)| {
                lc1.evaluate_hypercube(col1, bits1) * lc2.evaluate_hypercube(col2, bits2)
            })
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

    pub fn evaluate_hypercube(
        &self,
        col1: &[BaseField],
        bits1: usize,
        col2: &[BaseField],
        bits2: usize,
    ) -> BaseField {
        let random_k = &self.random_k;
        self.set
            .expressions
            .iter()
            .enumerate()
            .map(|(k, expression)| {
                let constraint_mask = Var::<K>(k).evaluate(random_k);
                constraint_mask * expression.evaluate_hypercube(col1, bits1, col2, bits2)
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

    pub fn evaluate_hypercube(
        &self,
        row: &[BaseField],
        row_bits: usize,
        col: &[BaseField],
        col_bits: usize,
    ) -> BaseField {
        let mut res = BaseField::from(0);
        for (i, coeffs) in self.rows.iter().enumerate() {
            let row_mask = Var::<I>(i).evaluate_hypercube(row, row_bits);
            for (j, coeff) in coeffs.iter().enumerate() {
                let col_mask = Var::<J>(j).evaluate_hypercube(col, col_bits);
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
        let one = BaseField::from(1);
        let pass = |i| {
            let a = self.data[i];
            let b = b[i];
            let c = c[i];
            a * b * c + (one - a) * (one - b) * (one - c)
        };
        (0..I).map(pass).product()
    }

    // Assumes `b.len() == c.len()`. No need for two bits, since they must be equal
    pub fn evaluate_hypercube(&self, b: &[BaseField], c: &[BaseField], bits: usize) -> BaseField {
        let one = BaseField::from(1);
        let pass = |i| {
            let a = self.data[i];
            let b = b[i];
            let c = c[i];
            a * b * c + (one - a) * (one - b) * (one - c)
        };
        let mut acc = (0..b.len()).map(pass).product();
        for i in 0..(I - b.len()) {
            let a = self.data[I - 1 - i];
            if (bits >> i) & 1 == 1 {
                acc *= a;
            } else {
                acc *= one - a;
            }
        }
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

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

    fn to_hypercube(index: u64, len: usize) -> Vec<BaseField> {
        let mut points = vec![BaseField::from(0); len];
        (0..len).for_each(|i| {
            if (index >> i) & 1 == 1 {
                points[len - 1 - i] = BaseField::from(1);
            }
        });
        points
    }

    #[test]
    fn var_hypercube_test() {
        let mut rng = rand::thread_rng();
        const J: usize = 16;
        const L: usize = 7;
        let var = Var::<J>::new(rng.gen_range(0..1 << J));
        let mask = (1 << (J - L)) - 1;
        let bits = var.0 & mask;
        let points = [(); L].map(|_| BaseField::from(rng.gen::<u64>()));
        let eval_hypercube = var.evaluate_hypercube(&points, bits);

        let mut points = points.to_vec();
        points.extend(to_hypercube(bits as u64, J - L));
        let points = (&points[..]).try_into().unwrap();
        let eval = var.evaluate(points);

        assert_eq!(eval_hypercube, eval);
    }

    #[test]
    fn delta_hypercube_test() {
        let mut rng = rand::thread_rng();
        const J: usize = 16;
        const L: usize = 7;
        let delta = Delta {
            data: [(); J].map(|_| BaseField::from(rng.gen::<u64>())),
        };
        let points1 = [(); L].map(|_| BaseField::from(rng.gen::<u64>()));
        let points2 = [(); L].map(|_| BaseField::from(rng.gen::<u64>()));
        let bits = rng.gen_range(0..1 << (J - L));
        let eval_hypercube = delta.evaluate_hypercube(&points1, &points2, bits);

        let mut points1 = points1.to_vec();
        points1.extend(to_hypercube(bits as u64, J - L));
        let points1 = (&points1[..]).try_into().unwrap();
        let mut points2 = points2.to_vec();
        points2.extend(to_hypercube(bits as u64, J - L));
        let points2 = (&points2[..]).try_into().unwrap();
        let eval = delta.evaluate(points1, points2);

        assert_eq!(eval_hypercube, eval);
    }
}
