use std::ops::{Add, Mul, Sub};

use crate::{
    constraints::{LinearCombination, Quadratic, Var},
    fields::BaseField,
};

pub enum Expr<const J: usize> {
    Elem(BaseField),
    Var(Var<J>),
    Add(Box<Expr<J>>, Box<Expr<J>>),
    Sub(Box<Expr<J>>, Box<Expr<J>>),
    Mul(Box<Expr<J>>, Box<Expr<J>>),
}

impl<const J: usize> Expr<J> {
    pub fn var(index: usize) -> Expr<J> {
        Expr::Var(Var::new(index))
    }

    pub fn to_linear_combination(&self) -> Option<LinearCombination<J>> {
        self.to_linear_combination_acc(BaseField::from(1))
    }

    fn to_linear_combination_acc(&self, acc: BaseField) -> Option<LinearCombination<J>> {
        match self {
            Expr::Var(var) => {
                let lc = LinearCombination([(acc, *var)].into());
                Some(lc)
            }
            Expr::Mul(a, b) => {
                if let &Expr::Elem(c) = a.as_ref() {
                    return b.to_linear_combination_acc(c * acc);
                }
                if let &Expr::Elem(c) = b.as_ref() {
                    return a.to_linear_combination_acc(c * acc);
                }
                None
            }
            Expr::Add(a, b) => {
                let a = a.to_linear_combination_acc(acc)?.0;
                let b = b.to_linear_combination_acc(acc)?.0;
                let mut c = Vec::with_capacity(a.len() + b.len());
                c.extend(a);
                c.extend(b);
                Some(LinearCombination(c.into()))
            }
            Expr::Sub(a, b) => {
                let a = a.to_linear_combination_acc(acc)?.0;
                let b = b.to_linear_combination_acc(-acc)?.0;
                let mut c = Vec::with_capacity(a.len() + b.len());
                c.extend(a);
                c.extend(b);
                Some(LinearCombination(c.into()))
            }
            Expr::Elem(..) => None,
        }
    }

    pub fn to_quadratic(&self) -> Option<Quadratic<J>> {
        self.to_quadratic_acc(BaseField::from(1))
    }

    fn to_quadratic_acc(&self, acc: BaseField) -> Option<Quadratic<J>> {
        match self {
            Expr::Add(a, b) => {
                let a = a.to_quadratic_acc(acc)?.0;
                let b = b.to_quadratic_acc(acc)?.0;
                let mut c = Vec::with_capacity(a.len() + b.len());
                c.extend(a);
                c.extend(b);
                Some(Quadratic(c.into()))
            }
            Expr::Sub(a, b) => {
                let a = a.to_quadratic_acc(acc)?.0;
                let b = b.to_quadratic_acc(-acc)?.0;
                let mut c = Vec::with_capacity(a.len() + b.len());
                c.extend(a);
                c.extend(b);
                Some(Quadratic(c.into()))
            }
            Expr::Mul(a, b) => {
                let a = a.to_linear_combination_acc(acc)?;
                let b = b.to_linear_combination()?;
                Some(Quadratic([(a, b)].into()))
            }
            Expr::Var(..) => None,
            Expr::Elem(..) => None,
        }
    }
}

impl<const J: usize> From<u64> for Expr<J> {
    fn from(elem: u64) -> Self {
        Expr::Elem(BaseField::from(elem))
    }
}

impl<const J: usize> Add for Expr<J> {
    type Output = Expr<J>;
    fn add(self, other: Expr<J>) -> Expr<J> {
        Expr::Add(self.into(), other.into())
    }
}

impl<const J: usize> Sub for Expr<J> {
    type Output = Expr<J>;
    fn sub(self, other: Expr<J>) -> Expr<J> {
        Expr::Sub(self.into(), other.into())
    }
}

impl<const J: usize> Mul for Expr<J> {
    type Output = Expr<J>;
    fn mul(self, other: Expr<J>) -> Expr<J> {
        Expr::Mul(self.into(), other.into())
    }
}
