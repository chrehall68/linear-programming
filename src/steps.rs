// convert a system of equations into canonical form
// for a given number of variables
// returns the intermediate steps too
use crate::problem::Tableau;
use fraction::{Fraction, Zero};

/// Convert a system of equations into canonical form
/// requires that the assigned basic variables are linearly independent
/// Returns [initial tableau, intermediate steps]. If the tableau is already in canonical form,
/// only returns [initial tableau]
fn canonicalize(mut prev_tableau: Tableau) -> Result<Vec<Tableau>, &'static str> {
    let mut steps = vec![prev_tableau.clone()];
    for pivot_row in 0..prev_tableau.constraints.len() {
        // find a row where the chosen variable is not zero
        // and we haven't already used that row as a pivot
        let tableau_row = &prev_tableau.constraints[pivot_row];
        let pivot_var = tableau_row.basic_var;
        let pivot_val = tableau_row.row.coefficients[pivot_var];
        if pivot_val == Fraction::zero() {
            return Err("Variables are not linearly independent");
        }
        let normalized_row = tableau_row.row.clone() / pivot_val;
        // now normalize the pivot row
        let mut next_tableau = prev_tableau.clone();
        next_tableau.constraints[pivot_row].row = normalized_row.clone();
        // and subtract it from the other rows
        for i in 0..prev_tableau.constraints.len() {
            if i == pivot_row {
                continue;
            }
            let scaler = next_tableau.constraints[i].row.coefficients[pivot_var];
            next_tableau.constraints[i].row -= normalized_row.clone() * scaler;
        }
        // and from the objective row
        let objective_scaler = next_tableau.objective.coefficients[pivot_var];
        next_tableau.objective -= normalized_row * objective_scaler;
        // and add it to the result and advance
        if prev_tableau != next_tableau {
            steps.push(next_tableau.clone());
        }
        prev_tableau = next_tableau;
    }
    Ok(steps)
}

/// Attempt to solve a linear program by using the simplex method
/// given an initial tableau.
/// The given variables in the initial tableau must form a basic feasible solution,
/// and if they don't, this will panic
pub fn simplex_method(initial_tableau: Tableau) -> Result<Vec<Tableau>, &'static str> {
    let mut steps = canonicalize(initial_tableau)?;
    let mut last = steps.last().unwrap().clone();

    // first make sure we're in a feasible solution
    // if we are, then all other things we visit will also be feasible
    // simplex method requires all RHS >= 0
    // and no restriction on whether objective coefficients are >= 0
    for table_row in &last.constraints {
        if table_row.row.rhs < Fraction::zero() {
            return Err("Cannot use simplex method; RHS must be >= 0");
        }
    }
    // then we just need to check the objective
    let find_most_negative_var = |last: &Tableau| -> Option<(usize, Fraction)> {
        last.objective
            .coefficients
            .iter()
            .enumerate()
            .min_by_key(|(_index, val)| *val)
            .map(|(idx, &val)| (idx, val))
    };
    while let Some((to_enter, _)) = find_most_negative_var(&last) {
        // make sure we're in a feasible solution
        // there's at least one variable with a negative coefficient in our objective
        // so now verify that the objective is bounded
        for (var, &coeff) in last.objective.coefficients.iter().enumerate() {
            if coeff < Fraction::zero() {
                // make sure that the variable is bounded
                let bounded = last
                    .constraints
                    .iter()
                    .any(|table_row| table_row.row.coefficients[var] > Fraction::zero());
                if !bounded {
                    return Err("Objective is unbounded");
                }
            }
        }
        // then see what variable needs to exit the bfs
        let to_exit_row = last
            .constraints
            .iter()
            .enumerate()
            .filter_map(|(row_idx, table_row)| {
                let row = &table_row.row;
                if row.coefficients[to_enter] > Fraction::zero() {
                    Some((row_idx, row.rhs / row.coefficients[to_enter]))
                } else {
                    None
                }
            })
            .min_by_key(|(_i, val)| *val)
            .unwrap()
            .0;
        last.constraints[to_exit_row].basic_var = to_enter;
        // canonicalize with respect to the new basic variable, and then advance
        let mut next_steps = canonicalize(last.clone())?;
        steps.append(&mut next_steps);
        last = steps.last().unwrap().clone();
    }
    Ok(steps)
}

/// Attempt to solve a linear program by using the simplex method
/// given an initial tableau
/// The given variables in the initial tableau must form a basic feasible solution,
/// and if they don't, this will panic
pub fn dual_simplex_method(initial_tableau: Tableau) -> Result<Vec<Tableau>, &'static str> {
    let mut steps = canonicalize(initial_tableau)?;
    let mut last = steps.last().unwrap().clone();

    // first make sure we're in a feasible solution
    // if we are, then all other things we visit will also be feasible
    // dual simplex requires that all objective coefficients are >= 0
    // and no restriction on whether RHS is >= 0
    for coeff in &last.objective.coefficients {
        if *coeff < Fraction::zero() {
            return Err("Cannot use dual simplex; objective coefficients must be >= 0");
        }
    }
    // then we just need to check the objective
    let find_most_negative_row = |last: &Tableau| -> Option<(usize, Fraction)> {
        last.constraints
            .iter()
            .enumerate()
            .min_by_key(|(_index, table_row)| table_row.row.rhs)
            .map(|(idx, table_row)| (idx, table_row.row.rhs))
    };
    while let Some((to_exit_row, _)) = find_most_negative_row(&last) {
        // make sure we're in a feasible solution
        for row in last.constraints.iter() {
            if row.row.rhs < Fraction::zero() {
                // make sure that the row is satisfiable (ie there exists a coefficient < 0)
                let bounded = row
                    .row
                    .coefficients
                    .iter()
                    .any(|&coeff| coeff < Fraction::zero());
                if !bounded {
                    return Err("No feasible solution");
                }
            }
        }

        // then see what variable needs to enter the bfs
        let to_enter = last.constraints[to_exit_row]
            .row
            .coefficients
            .iter()
            .enumerate()
            .filter_map(|(row_idx, &coeff)| {
                if coeff < Fraction::zero() {
                    Some((row_idx, last.objective.coefficients[row_idx] / coeff))
                } else {
                    None
                }
            })
            .max_by_key(|(_i, val)| *val)
            .unwrap()
            .0;
        last.constraints[to_exit_row].basic_var = to_enter;
        // canonicalize with respect to the new basic variable, and then advance
        let mut next_steps = canonicalize(last.clone())?;
        steps.append(&mut next_steps);
        last = steps.last().unwrap().clone();
    }
    Ok(steps)
}
