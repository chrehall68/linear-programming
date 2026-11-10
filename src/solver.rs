use fraction::{Fraction, One, Zero};

use crate::{
    problem::{Inequality, Problem, Row, TableRow, Tableau},
    steps::{canonicalize, simplex_method},
};

/// This function will attempt to solve a linear program
/// by running the two-pass simplex method:
/// The first pass finds a BFS using artificial variables, and the second pass finds a solution
/// using the BFS from the first pass.
///
/// If there is no solution or the objective is unbounded, it will error
/// All RHS values must be >= 0
///
/// Returns (first_pass_steps, second_pass_steps)
pub fn solve(mut problem: Problem) -> Result<(Vec<Tableau>, Vec<Tableau>), &'static str> {
    // first add the slack and artificial variables
    let initial_vars = problem.num_vars();
    let mut added_slacks = 0;
    let mut artificial_idxs = Vec::new();
    let mut variable_names: Vec<_> = (0..initial_vars).map(|v| format!("x{}", v + 1)).collect();
    let mut tableau_rows: Vec<_> = problem
        .constraints
        .iter_mut()
        .map(|constraint| {
            // first swap sign if necessary
            if constraint.row.rhs < Fraction::zero() {
                constraint.row *= -Fraction::one();
                constraint.inequality = constraint.inequality.opposite();
            }
            // add extras
            for _ in 0..(added_slacks + artificial_idxs.len()) {
                constraint.row.coefficients.push(Fraction::zero());
            }
            // then maybe add slack or artificial or both
            // and have whatever we use be the last variable in the row
            match constraint.inequality {
                Inequality::Le => {
                    // since we have ... <= b, just do
                    // ... + s = b
                    constraint.row.coefficients.push(Fraction::one());
                    variable_names.push(format!("s{}", added_slacks + 1));
                    added_slacks += 1;
                }
                Inequality::Eq => {
                    // since we have ... = b, do
                    // ... + a = b, and we'd expect a to become 0
                    constraint.row.coefficients.push(Fraction::one());
                    variable_names.push(format!("a{}", artificial_idxs.len() + 1));
                    artificial_idxs.push(constraint.row.coefficients.len() - 1);
                }
                Inequality::Ge => {
                    // we have ... >= b, so we make it
                    // ... - s = b, but then we also insert an artificial to use
                    // as a pivot, giving
                    // ... - s + a = b
                    // slack first
                    constraint.row.coefficients.push(-Fraction::one());
                    variable_names.push(format!("s{}", added_slacks + 1));
                    added_slacks += 1;
                    // then artificial since that's what we'll use as our BV
                    constraint.row.coefficients.push(Fraction::one());
                    artificial_idxs.push(constraint.row.coefficients.len() - 1);
                    variable_names.push(format!("a{}", artificial_idxs.len() + 1));
                }
            }
            TableRow {
                row: constraint.row.clone(),
                basic_var: constraint.row.coefficients.len() - 1,
            }
        })
        .collect();
    // pad remaining
    for tableau_row in tableau_rows.iter_mut() {
        while tableau_row.row.coefficients.len()
            < initial_vars + added_slacks + artificial_idxs.len()
        {
            tableau_row.row.coefficients.push(Fraction::zero());
        }
    }
    // then construct the objective row
    let mut objective_coeffs =
        vec![Fraction::zero(); initial_vars + added_slacks + artificial_idxs.len()];
    for &artificial_var in &artificial_idxs {
        objective_coeffs[artificial_var] = Fraction::one();
    }
    let artificial_tableau = Tableau {
        objective: Row {
            coefficients: objective_coeffs,
            rhs: Fraction::zero(),
        },
        constraints: tableau_rows,
    };

    // then run simplex method
    let mut artificial_steps = simplex_method(artificial_tableau)?;
    // we know that it's possible iff the objective was minimized to 0
    if artificial_steps.last().unwrap().objective.rhs != Fraction::zero() {
        return Err("No feasible solution");
    }

    // now construct tableau w/o artificials
    // if there's any remaining artificials, then we can
    // pivot on a nonzero regular variable to move the artificial out of our BFS
    // and if there are no nonzero regular variables in the row, the row is redundant and can be skipped
    let mut last = artificial_steps.last().unwrap().clone();
    let mut is_artificial = vec![false; initial_vars + added_slacks + artificial_idxs.len()];
    let mut removed_row = vec![false; last.constraints.len()];
    for &artificial_var in &artificial_idxs {
        is_artificial[artificial_var] = true;
    }
    while let Some(row_idx) = last
        .constraints
        .iter()
        .enumerate()
        .filter_map(|(row_idx, tableau_row)| {
            if !removed_row[row_idx] && !is_artificial[tableau_row.basic_var] {
                Some(row_idx)
            } else {
                None
            }
        })
        .next()
    {
        let pivot_var = (0..(initial_vars + added_slacks + artificial_idxs.len())).find(|&var| {
            !is_artificial[var]
                && last.constraints[row_idx].row.coefficients[var] != Fraction::zero()
        });
        if let Some(pivot_var) = pivot_var {
            // pivot on this
            last.constraints[row_idx].basic_var = pivot_var;
            let mut pivot_steps = canonicalize(last)?;
            artificial_steps.append(&mut pivot_steps);
            last = artificial_steps.last().unwrap().clone();
        } else {
            removed_row[row_idx] = true;
        }
    }

    // remove redundant rows
    let filtered_constraints: Vec<_> = last
        .constraints
        .iter()
        .enumerate()
        .filter_map(|(row_idx, tableau_row)| {
            if !removed_row[row_idx] {
                Some(tableau_row.clone())
            } else {
                None
            }
        })
        .collect();

    // now solve the problem with the reduced tableau and original objective
    let mut original_objective = problem.objective.clone();
    while original_objective.coefficients.len()
        < initial_vars + added_slacks + artificial_idxs.len()
    {
        original_objective.coefficients.push(Fraction::zero());
    }

    // before starting, be careful about 0 constraints
    if filtered_constraints.len() == 0 {
        if original_objective
            .coefficients
            .iter()
            .all(|x| *x > Fraction::zero())
        {
            // one solution: have all variables = 0
            return Ok((
                artificial_steps,
                vec![Tableau {
                    objective: Row {
                        coefficients: vec![
                            Fraction::zero();
                            initial_vars + added_slacks + artificial_idxs.len()
                        ],
                        rhs: Fraction::zero(),
                    },
                    constraints: vec![TableRow {
                        row: Row {
                            coefficients: vec![
                                Fraction::zero();
                                initial_vars + added_slacks + artificial_idxs.len()
                            ],
                            rhs: Fraction::zero(),
                        },
                        basic_var: 0,
                    }],
                }],
            ));
        } else {
            // infinitely many solutions since there's at least one variable
            // with objective value <= 0. Then, we can just take any nonnegative real value for that
            // to get infinite solutions
            return Err(
                "Unbounded since there are no constraints and at least one variable with objective value <= 0",
            );
        }
    }

    // otherwise, we have at least one constraint
    // so we can solve it regularly
    let second_pass_tableau = Tableau {
        objective: original_objective,
        constraints: filtered_constraints,
    };
    let second_pass_steps = simplex_method(second_pass_tableau)?;
    Ok((artificial_steps, second_pass_steps))
}
