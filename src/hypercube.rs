use crate::{
    constraints::{ConstraintSetMatrix, Delta, LinearCombination, Quadratic, Trace, Var, BQCS},
    fields::BaseField,
};
use std::array::from_fn;

// All functions assume `bits` is within the multilinear's range, and that `I`,`J` < `usize::BITS`
impl<const J: usize> Var<J> {
    #[inline]
    pub fn evaluate_hypercube(&self, bits: usize) -> BaseField {
        if bits == self.0 {
            return BaseField::from(1);
        }
        BaseField::from(0)
    }
}

impl<const J: usize> LinearCombination<J> {
    pub fn evaluate_hypercube(&self, bits: usize) -> BaseField {
        self.0
            .iter()
            .map(|(c, x)| *c * x.evaluate_hypercube(bits))
            .sum()
    }
}

impl<const J: usize> Quadratic<J> {
    pub fn evaluate_hypercube(&self, bits1: usize, bits2: usize) -> BaseField {
        self.0
            .iter()
            .map(|(lc1, lc2)| lc1.evaluate_hypercube(bits1) * lc2.evaluate_hypercube(bits2))
            .sum()
    }
}

impl<const J: usize, const K: usize> ConstraintSetMatrix<J, K>
where
    [(); 1 << K]:,
{
    pub fn evaluate_hypercube(&self, bits1: usize, bits2: usize) -> BaseField {
        let constraint_masks = self.constraint_masks();
        self.set
            .expressions
            .iter()
            .zip(constraint_masks)
            .map(|(expression, &constraint_mask)| {
                constraint_mask * expression.evaluate_hypercube(bits1, bits2)
            })
            .sum()
    }

    // TODO: Build tables for the linear combinations instead, so that it's not `J^2`
    pub fn build_table(&self) -> Box<[[BaseField; 1 << J]; 1 << J]> {
        Box::new(from_fn(|j1| from_fn(|j2| self.evaluate_hypercube(j1, j2))))
    }
}

impl<const I: usize, const J: usize> Trace<I, J>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
{
    pub fn evaluate_hypercube(&self, row_bits: usize, col_bits: usize) -> BaseField {
        self.rows[row_bits][col_bits]
    }

    pub fn build_table(&self) -> Box<[[BaseField; 1 << J]; 1 << I]> {
        self.rows.clone()
    }
}

impl<const I: usize> Delta<I> {
    // No need for two bits, since they must be equal
    pub fn evaluate_hypercube(&self, bits: usize) -> BaseField {
        let one = BaseField::from(1);
        let mut acc = BaseField::from(1);
        for i in 0..I {
            let a = self.data[I - 1 - i];
            if (bits >> i) & 1 == 1 {
                acc *= a;
            } else {
                acc *= one - a;
            }
        }
        acc
    }

    pub fn build_table(&self) -> Box<[BaseField; 1 << I]> {
        let vec: Vec<_> = (0..1 << I).map(|i| self.evaluate_hypercube(i)).collect();
        vec.try_into().unwrap()
    }
}

impl<const I: usize, const J: usize, const K: usize> BQCS<I, J, K>
where
    [(); 1 << I]:,
    [(); 1 << J]:,
    [(); 1 << K]:,
{
    pub fn evaluate_hypercube(
        &self,
        row_bits: usize,
        col1_bits: usize,
        col2_bits: usize,
    ) -> BaseField {
        let d = self.delta.evaluate_hypercube(row_bits);
        let a = self.matrix.evaluate_hypercube(col1_bits, col2_bits);
        let w1 = self.trace.evaluate_hypercube(row_bits, col1_bits);
        let w2 = self.trace.evaluate_hypercube(row_bits, col2_bits);
        d * a * w1 * w2
    }
}
