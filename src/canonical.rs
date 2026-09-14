// convert a system of equations into canonical form
// for a given number of variables
// returns the intermediate steps too
use fraction::{Fraction, One, Zero};

pub struct Step {
    /// (coefficients, RHS)
    pub equations: Vec<(Vec<Fraction>, Fraction)>,
    /// (coefficients, RHS)
    pub objective: (Vec<Fraction>, Fraction),
    /// which row was used as a pivot, 0-indexed
    pub pivot_row: usize, // assumes that objective is the last equation
                          // and that the rhs is the last coefficient
}
impl Step {
    fn extract(mut coefficients: Vec<Fraction>) -> (Vec<Fraction>, Fraction) {
        let rhs = coefficients.pop().unwrap();
        (coefficients, rhs)
    }
    fn from_equations(mut coefficients: Vec<Vec<Fraction>>, pivot_row: usize) -> Self {
        let objective = Self::extract(coefficients.pop().unwrap());
        let equations: Vec<_> = coefficients
            .into_iter()
            .map(|coefficients| Self::extract(coefficients))
            .collect();
        Step {
            equations,
            objective,
            pivot_row,
        }
    }
}
pub fn canonicalize(
    equations: Vec<(Vec<Fraction>, Fraction)>,
    mut objective: Vec<Fraction>,
    variables: Vec<usize>,
) -> Vec<Step> {
    let num_vars = equations[0].0.len();
    let num_equations = equations.len();
    assert!(num_equations == variables.len());
    assert!(num_vars >= variables.len());
    let mut was_used = vec![false; num_vars];
    let mut result = Vec::new();
    // store everything just as a vector of fractions for ease of use
    let mut cur_equations: Vec<_> = equations
        .into_iter()
        .map(|(mut coeffs, rhs)| {
            coeffs.push(rhs);
            coeffs
        })
        .collect();
    objective.push(Fraction::zero());
    cur_equations.push(objective);
    // then just use each of the variables
    for chosen in variables {
        // find a row where the chosen variable is not zero
        // and we haven't already used that row as a pivot
        let mut pivot_row = None;
        for i in 0..num_equations {
            // don't need to check the objective
            if cur_equations[i][chosen] != Fraction::zero() && !was_used[i] {
                pivot_row = Some(i);
                break;
            }
        }
        assert!(pivot_row.is_some());
        was_used[pivot_row.unwrap()] = true;
        let pivot_row = pivot_row.unwrap();
        let pivot_val = cur_equations[pivot_row][chosen];
        // now normalize the pivot row
        let pivot_inv = Fraction::one() / pivot_val;
        let mut next_equations = cur_equations.clone();
        next_equations[pivot_row] = next_equations[pivot_row]
            .iter()
            .map(|x| x * pivot_inv)
            .collect();
        // and subtract it from the other rows
        for i in 0..(num_equations + 1) {
            if i == pivot_row {
                continue;
            }
            let scaler = next_equations[i][chosen];
            next_equations[i] = next_equations[i]
                .iter()
                .zip(next_equations[pivot_row].iter())
                .map(|(cur_val, pivot_val)| cur_val - scaler * *pivot_val)
                .collect();
        }
        // and add it to the result and advance
        result.push(Step::from_equations(next_equations.clone(), pivot_row));
        cur_equations = next_equations;
    }
    result
}
