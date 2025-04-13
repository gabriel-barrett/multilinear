use ark_ff::AdditiveGroup;

use crate::{
    constraints::{ConstraintSetMatrix, Delta, LinearCombination, Quadratic, Trace, Var, BQCS},
    fields::BaseField,
};

impl<const J: usize> Var<J> {
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

impl<const J: usize> LinearCombination<J> {
    pub fn evaluate_hypercube(&self, points: &[BaseField], bits: usize) -> BaseField {
        self.0
            .iter()
            .map(|(c, x)| *c * x.evaluate_hypercube(points, bits))
            .sum()
    }
}

impl<const J: usize> Quadratic<J> {
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

impl<const J: usize, const K: usize> ConstraintSetMatrix<J, K>
where
    [(); 1 << K]:,
{
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

impl<const I: usize, const J: usize> Trace<I, J>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    pub fn evaluate_hypercube(
        &self,
        row: &[BaseField],
        row_bits: usize,
        col: &[BaseField],
        col_bits: usize,
    ) -> BaseField {
        let mut res = BaseField::from(0);
        let row_jump = 1 << (I - row.len());
        let col_jump = 1 << (J - col.len());
        for i in 0..1 << row.len() {
            let i = i * row_jump + row_bits;
            let coeffs = &self.rows[i];
            let row_mask = Var::<I>(i).evaluate_hypercube(row, row_bits);
            for j in 0..1 << col.len() {
                let j = j * col_jump + col_bits;
                let col_mask = Var::<J>(j).evaluate_hypercube(col, col_bits);
                res += coeffs[j] * row_mask * col_mask;
            }
        }
        res
    }

    pub fn evaluate_hypercube_row(
        &self,
        row: usize,
        col: &[BaseField],
        col_bits: usize,
    ) -> BaseField {
        let mut res = BaseField::from(0);
        let coeffs = self.rows[row];
        for (j, coeff) in coeffs.iter().enumerate() {
            let col_mask = Var::<J>(j).evaluate_hypercube(col, col_bits);
            res += *coeff * col_mask;
        }
        res
    }

    pub fn evaluate_hypercube_col(
        &self,
        row: &[BaseField],
        row_bits: usize,
        col: usize,
    ) -> BaseField {
        let mut res = BaseField::from(0);
        for (i, coeffs) in self.rows.iter().enumerate() {
            let row_mask = Var::<I>(i).evaluate_hypercube(row, row_bits);
            res += coeffs[col] * row_mask;
        }
        res
    }
}

impl<const I: usize> Delta<I> {
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

impl<const I: usize, const J: usize, const K: usize> BQCS<I, J, K>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
    [(); 1 << K]:,
{
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_hypercube(
        &self,
        row1: &[BaseField],
        row2: &[BaseField],
        row_bits: usize,
        col1: &[BaseField],
        col1_bits: usize,
        col2: &[BaseField],
        col2_bits: usize,
    ) -> BaseField {
        let d = self.delta.evaluate_hypercube(row1, row2, row_bits);
        let a = self
            .matrix
            .evaluate_hypercube(col1, col1_bits, col2, col2_bits);
        let c = d * a;
        if c == BaseField::ZERO {
            return BaseField::ZERO;
        }
        let w1 = self
            .trace
            .evaluate_hypercube(row1, row_bits, col1, col1_bits);
        let w2 = self
            .trace
            .evaluate_hypercube(row2, row_bits, col2, col2_bits);
        c * w1 * w2
    }

    pub fn evaluate_hypercube_row(
        &self,
        row: usize,
        col1: &[BaseField],
        col1_bits: usize,
        col2: &[BaseField],
        col2_bits: usize,
    ) -> BaseField {
        let d = self.delta.evaluate_hypercube(&[], &[], row);
        let a = self
            .matrix
            .evaluate_hypercube(col1, col1_bits, col2, col2_bits);
        let c = d * a;
        if c == BaseField::ZERO {
            return BaseField::ZERO;
        }
        let w1 = self.trace.evaluate_hypercube_row(row, col1, col1_bits);
        let w2 = self.trace.evaluate_hypercube_row(row, col2, col2_bits);
        c * w1 * w2
    }

    pub fn evaluate_hypercube_col(
        &self,
        row1: &[BaseField],
        row2: &[BaseField],
        row_bits: usize,
        col1: usize,
        col2: usize,
    ) -> BaseField {
        let d = self.delta.evaluate_hypercube(row1, row2, row_bits);
        let a = self.matrix.evaluate_hypercube(&[], col1, &[], col2);
        let c = d * a;
        if c == BaseField::ZERO {
            return BaseField::ZERO;
        }
        let w1 = self.trace.evaluate_hypercube_col(row1, row_bits, col1);
        let w2 = self.trace.evaluate_hypercube_col(row2, row_bits, col2);
        c * w1 * w2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::random_base;
    use rand::Rng;

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
        let rng = &mut rand::rng();
        const J: usize = 16;
        const L: usize = 7;
        let var = Var::<J>::new(rng.random_range(0..1 << J));
        let mask = (1 << (J - L)) - 1;
        let bits = var.0 & mask;
        let points = [(); L].map(|_| random_base(rng));
        let eval_hypercube = var.evaluate_hypercube(&points, bits);

        let mut points = points.to_vec();
        points.extend(to_hypercube(bits as u64, J - L));
        let points = (&points[..]).try_into().unwrap();
        let eval = var.evaluate(points);

        assert_eq!(eval_hypercube, eval);
    }

    #[test]
    fn delta_hypercube_test() {
        let rng = &mut rand::rng();
        const J: usize = 16;
        const L: usize = 7;
        let delta = Delta {
            data: [(); J].map(|_| random_base(rng)),
        };
        let points1 = [(); L].map(|_| random_base(rng));
        let points2 = [(); L].map(|_| random_base(rng));
        let bits = rng.random_range(0..1 << (J - L));
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
