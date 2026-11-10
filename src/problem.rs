use fraction::Fraction;
use std::ops::{AddAssign, Div, DivAssign, Mul, MulAssign, SubAssign};

// ==================================================
// Row Setup / Operations
// ==================================================
#[derive(Clone, PartialEq, Eq)]
pub struct Row {
    pub coefficients: Vec<Fraction>,
    pub rhs: Fraction,
}
impl AddAssign for Row {
    fn add_assign(&mut self, rhs: Row) {
        assert_eq!(self.coefficients.len(), rhs.coefficients.len());
        self.coefficients
            .iter_mut()
            .zip(&rhs.coefficients)
            .for_each(|(a, b)| *a += b);
        self.rhs += rhs.rhs;
    }
}
impl SubAssign for Row {
    fn sub_assign(&mut self, rhs: Row) {
        assert_eq!(self.coefficients.len(), rhs.coefficients.len());
        self.coefficients
            .iter_mut()
            .zip(&rhs.coefficients)
            .for_each(|(a, b)| *a -= b);
        self.rhs -= rhs.rhs;
    }
}
impl DivAssign<Fraction> for Row {
    fn div_assign(&mut self, rhs: Fraction) {
        self.coefficients.iter_mut().for_each(|c| *c /= rhs);
        self.rhs /= rhs;
    }
}
impl MulAssign<Fraction> for Row {
    fn mul_assign(&mut self, rhs: Fraction) {
        self.coefficients.iter_mut().for_each(|c| *c *= rhs);
        self.rhs *= rhs;
    }
}
impl Mul<Fraction> for Row {
    type Output = Row;
    fn mul(self, rhs: Fraction) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}
impl Div<Fraction> for Row {
    type Output = Row;
    fn div(self, rhs: Fraction) -> Self::Output {
        let mut output = self.clone();
        output /= rhs;
        output
    }
}
// ==================================================
// Problem Setup
// ==================================================
pub enum Inequality {
    Eq,
    Le,
    Ge,
}
impl Inequality {
    pub fn opposite(&self) -> Inequality {
        match *self {
            Inequality::Eq => Inequality::Eq,
            Inequality::Le => Inequality::Ge,
            Inequality::Ge => Inequality::Le,
        }
    }
}
pub struct InputConstraint {
    pub row: Row,
    pub inequality: Inequality,
}
pub struct Problem {
    pub constraints: Vec<InputConstraint>,
    pub objective: Row,
}
impl Problem {
    pub fn new(constraints: Vec<InputConstraint>, objective: Row) -> Problem {
        // check that we have at least one constraint
        assert!(!constraints.is_empty(), "must have at least one constraint");
        // check that all constraints have the same number of coefficients
        let num_vars = constraints[0].row.coefficients.len();
        assert!(
            constraints[0].row.coefficients.len() > 0,
            "must have at least one variable"
        );
        for constraint in &constraints {
            assert_eq!(
                constraint.row.coefficients.len(),
                num_vars,
                "all constraints must have same number of coefficients"
            );
        }
        assert_eq!(
            objective.coefficients.len(),
            num_vars,
            "objective must have same number of coefficients as constraints"
        );
        Problem {
            constraints,
            objective,
        }
    }
    pub(crate) fn num_vars(&self) -> usize {
        self.constraints[0].row.coefficients.len()
    }
    pub(crate) fn num_constraints(&self) -> usize {
        self.constraints.len()
    }
}

// ==================================================
// Answer Tableau
// ==================================================
/// A row in the answer tableau
/// Like a regular Row, but with an associated basic variable
#[derive(Clone, PartialEq, Eq)]
pub struct TableRow {
    pub row: Row,
    /// 0-indexed index of the basic variable associated with this row
    pub basic_var: usize,
}
/// A single answer tableau
#[derive(Clone, PartialEq, Eq)]
pub struct Tableau {
    pub constraints: Vec<TableRow>,
    pub objective: Row,
}
