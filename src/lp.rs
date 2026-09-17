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
    /// which variable was used as a pivot
    pub pivot_var: usize,
}
impl Step {
    fn extract(mut coefficients: Vec<Fraction>) -> (Vec<Fraction>, Fraction) {
        let rhs = coefficients.pop().unwrap();
        (coefficients, rhs)
    }
    fn from_equations(
        mut coefficients: Vec<Vec<Fraction>>,
        pivot_row: usize,
        pivot_var: usize,
    ) -> Self {
        let objective = Self::extract(coefficients.pop().unwrap());
        let equations: Vec<_> = coefficients
            .into_iter()
            .map(|coefficients| Self::extract(coefficients))
            .collect();
        Step {
            equations,
            objective,
            pivot_row,
            pivot_var,
        }
    }
}
// returns (steps, v such that v[row] = which basic variable is in that row)
fn canonicalize_inner(
    mut equations: Vec<Vec<Fraction>>,
    variables: Vec<usize>,
) -> (Vec<Step>, Vec<usize>) {
    let mut steps = Vec::new();
    let mut which_vars = vec![None; equations.len() - 1];

    for chosen in variables {
        // find a row where the chosen variable is not zero
        // and we haven't already used that row as a pivot
        let mut pivot_row = None;
        for i in 0..(equations.len() - 1) {
            // don't need to check the objective
            if equations[i][chosen] != Fraction::zero() && which_vars[i].is_none() {
                pivot_row = Some(i);
                break;
            }
        }
        assert!(pivot_row.is_some());
        which_vars[pivot_row.unwrap()] = Some(chosen);
        let pivot_row = pivot_row.unwrap();
        let pivot_val = equations[pivot_row][chosen];
        // now normalize the pivot row
        let pivot_inv = Fraction::one() / pivot_val;
        let mut next_equations = equations.clone();
        next_equations[pivot_row] = next_equations[pivot_row]
            .iter()
            .map(|x| x * pivot_inv)
            .collect();
        // and subtract it from the other rows
        for i in 0..equations.len() {
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
        steps.push(Step::from_equations(
            next_equations.clone(),
            pivot_row,
            chosen,
        ));
        equations = next_equations;
    }
    (steps, which_vars.into_iter().map(|x| x.unwrap()).collect())
}
fn preprocess(
    equations: Vec<(Vec<Fraction>, Fraction)>,
    mut objective: Vec<Fraction>,
    objective_val: Fraction,
    chosen_vars: usize,
) -> Vec<Vec<Fraction>> {
    let num_vars = equations[0].0.len();
    let num_equations = equations.len();
    assert!(num_equations == chosen_vars);
    assert!(num_vars >= chosen_vars);
    // store everything just as a vector of fractions for ease of use
    let mut cur_equations: Vec<_> = equations
        .into_iter()
        .map(|(mut coeffs, rhs)| {
            coeffs.push(rhs);
            coeffs
        })
        .collect();
    objective.push(objective_val);
    cur_equations.push(objective);
    cur_equations
}
pub fn canonicalize(
    equations: Vec<(Vec<Fraction>, Fraction)>,
    objective: Vec<Fraction>,
    variables: Vec<usize>,
) -> Vec<Step> {
    // then just use each of the variables
    let cur_equations = preprocess(equations, objective, Fraction::zero(), variables.len());
    canonicalize_inner(cur_equations, variables).0
}
/// attempt to solve a linear program
/// given an equation in standard form and a set of variables
/// to start with.
/// The given variables must form a basic feasible solution,
/// and if they don't, this will panic
pub fn solve(
    equations: Vec<(Vec<Fraction>, Fraction)>,
    objective: Vec<Fraction>,
    mut variables: Vec<usize>,
) -> Vec<Step> {
    let num_equations = equations.len();
    let cur_equations = preprocess(equations, objective, Fraction::zero(), variables.len());
    let (mut steps, mut which_vars) = canonicalize_inner(cur_equations, variables.clone());
    let mut last = steps.last().unwrap();
    let mut coeffs = last.equations.clone();
    let mut objective = last.objective.0.clone();

    // first make sure we're in a feasible solution
    // if we are, then all other things we visit will also be feasible
    for row in 0..num_equations {
        if coeffs[row].1 < Fraction::zero() {
            panic!("not in a feasible solution");
        }
    }
    // then we just need to check the objective
    while objective.iter().any(|x| x < &Fraction::zero()) {
        // make sure we're in a feasible solution
        // there's at least one variable with a negative coefficient in our objective
        // so now verify that the objective is bounded
        for v in 0..objective.len() {
            if objective[v] < Fraction::zero() {
                // make sure that the variable is bounded
                let bounded = (0..num_equations).any(|row| coeffs[row].0[v] > Fraction::zero());
                if !bounded {
                    panic!("objective is unbounded");
                }
            }
        }
        // objective is bounded, so now just find the most negative coefficient
        let most_negative = objective
            .iter()
            .enumerate()
            .min_by_key(|(_index, val)| *val)
            .unwrap()
            .0;
        // then see what variable gets replaced
        let to_replace_row = (0..num_equations)
            .filter_map(|row| {
                if coeffs[row].0[most_negative] > Fraction::zero() {
                    Some((row, coeffs[row].1 / coeffs[row].0[most_negative]))
                } else {
                    None
                }
            })
            .min_by_key(|(_i, val)| *val)
            .unwrap()
            .0;
        let to_replace_var = which_vars[to_replace_row];
        let mut next_variables: Vec<_> = variables
            .into_iter()
            .filter(|x| *x != to_replace_var)
            .collect();
        next_variables.push(most_negative);
        // run and then advance
        let next_equations = preprocess(coeffs, objective, last.objective.1, next_variables.len());
        let (mut next_steps, next_which_vars) =
            canonicalize_inner(next_equations, next_variables.clone());
        // only take the last one since all other variables are already
        // in canonical form, meaning the operations done on them are no-ops
        // and the only actual operation is the last operation, where we're pivoting
        // on the new variable
        steps.push(next_steps.pop().unwrap());
        which_vars = next_which_vars;
        variables = next_variables;
        last = steps.last().unwrap();
        coeffs = last.equations.clone();
        objective = last.objective.0.clone();
    }
    steps
}
// TODO - 2-phase simplex to find the BFS first
