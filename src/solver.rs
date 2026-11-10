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
/// Returns (first_pass_steps, second_pass_steps, variable_names)
pub fn solve(
    mut problem: Problem,
) -> Result<(Vec<Tableau>, Vec<Tableau>, Vec<String>), &'static str> {
    // first add the slack and artificial variables
    let initial_vars = problem.num_vars();
    let mut added_slacks = 0;
    let mut added_artificials = 0;
    // count initially that way we can construct all artificials at end
    let total_slacks = problem
        .constraints
        .iter()
        .filter(|c| c.inequality == Inequality::Le || c.inequality == Inequality::Ge)
        .count();
    let total_artificial = problem
        .constraints
        .iter()
        .filter(|c| c.inequality == Inequality::Eq || c.inequality == Inequality::Ge)
        .count();
    let mut variable_names: Vec<_> = (0..initial_vars).map(|v| format!("x{}", v + 1)).collect();
    for i in 1..=total_slacks {
        variable_names.push(format!("s{i}"));
    }
    for i in 1..=total_artificial {
        variable_names.push(format!("a{i}"));
    }
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
            for _ in 0..(total_slacks + total_artificial) {
                constraint.row.coefficients.push(Fraction::zero());
            }
            // then maybe add slack or artificial or both
            // and have whatever we use be the last variable in the row
            match constraint.inequality {
                Inequality::Le => {
                    // since we have ... <= b, just do
                    // ... + s = b
                    constraint.row.coefficients[initial_vars + added_slacks] = Fraction::one();
                    added_slacks += 1;
                }
                Inequality::Eq => {
                    // since we have ... = b, do
                    // ... + a = b, and we'd expect a to become 0
                    constraint.row.coefficients[initial_vars + total_slacks + added_artificials] =
                        Fraction::one();
                    added_artificials += 1;
                }
                Inequality::Ge => {
                    // we have ... >= b, so we make it
                    // ... - s = b, but then we also insert an artificial to use
                    // as a pivot, giving
                    // ... - s + a = b
                    // slack first
                    constraint.row.coefficients[initial_vars + added_slacks] = Fraction::one();
                    added_slacks += 1;
                    // then artificial since that's what we'll use as our BV
                    constraint.row.coefficients[initial_vars + total_slacks + added_artificials] =
                        Fraction::one();
                    added_artificials += 1;
                }
            }
            TableRow {
                row: constraint.row.clone(),
                basic_var: constraint.row.coefficients.len() - 1,
            }
        })
        .collect();
    // then construct the objective row
    let mut objective_coeffs =
        vec![Fraction::zero(); initial_vars + total_slacks + total_artificial];
    for i in 0..total_artificial {
        objective_coeffs[initial_vars + total_slacks + i] = Fraction::one();
    }
    let mut original_objective = problem.objective.clone();
    while original_objective.coefficients.len() < initial_vars + total_slacks + total_artificial {
        original_objective.coefficients.push(Fraction::zero());
    }
    let artificial_tableau = Tableau {
        objective: Row {
            coefficients: objective_coeffs,
            rhs: Fraction::zero(),
        },
        constraints: tableau_rows,
        secondary_objective: Some(original_objective),
    };

    // then run simplex method
    let mut artificial_steps = simplex_method(artificial_tableau)?;
    // we know that it's possible iff the objective was minimized to 0
    if artificial_steps.last().unwrap().objective.rhs != Fraction::zero() {
        return Err("No feasible solution");
    }

    // now construct tableau w/o artificials
    // if there's any remaining artificials in our BFS, then we can
    // pivot on a nonzero regular variable to move the artificial out of our BFS
    // and if there are no nonzero regular variables in the row, the row is redundant and can be skipped
    let mut last = artificial_steps.last().unwrap().clone();
    let mut removed_row = vec![false; last.constraints.len()];
    while let Some(row_idx) = last
        .constraints
        .iter()
        .enumerate()
        .filter_map(|(row_idx, tableau_row)| {
            if !removed_row[row_idx] && tableau_row.basic_var >= initial_vars + added_slacks {
                Some(row_idx)
            } else {
                None
            }
        })
        .next()
    {
        let pivot_var = (0..(initial_vars + added_slacks))
            .find(|&var| last.constraints[row_idx].row.coefficients[var] != Fraction::zero());
        if let Some(pivot_var) = pivot_var {
            // pivot on this
            last.constraints[row_idx].basic_var = pivot_var;
            let mut pivot_steps = canonicalize(last)?;
            artificial_steps.append(&mut pivot_steps);
            last = artificial_steps.last().unwrap().clone();
        } else {
            // safe to remove since we know rhs must be zero since total objective is 0
            removed_row[row_idx] = true;
        }
    }

    // remove redundant rows and truncate to just original variables and slacks
    let filtered_constraints: Vec<_> = last
        .constraints
        .iter()
        .enumerate()
        .filter_map(|(row_idx, tableau_row)| {
            if !removed_row[row_idx] {
                let mut tableau_row = tableau_row.clone();
                tableau_row
                    .row
                    .coefficients
                    .truncate(initial_vars + added_slacks);
                Some(tableau_row)
            } else {
                None
            }
        })
        .collect();

    // now solve the problem with the reduced tableau and original objective
    let mut next_objective = last.secondary_objective.unwrap();
    next_objective
        .coefficients
        .truncate(initial_vars + added_slacks);

    // before starting, be careful about 0 constraints
    if filtered_constraints.len() == 0 {
        if next_objective
            .coefficients
            .iter()
            .all(|x| *x >= Fraction::zero())
        {
            // one solution: have all variables = 0
            return Ok((
                artificial_steps,
                vec![Tableau {
                    objective: next_objective,
                    constraints: vec![TableRow {
                        row: Row {
                            coefficients: vec![Fraction::zero(); initial_vars + total_slacks],
                            rhs: Fraction::zero(),
                        },
                        basic_var: 0,
                    }],
                    secondary_objective: None,
                }],
                variable_names,
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
        objective: next_objective,
        constraints: filtered_constraints,
        secondary_objective: None,
    };
    let second_pass_steps = simplex_method(second_pass_tableau)?;
    Ok((artificial_steps, second_pass_steps, variable_names))
}
