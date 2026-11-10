use fraction::{Fraction, One, Zero};

use crate::{
    problem::{Inequality, Problem, Row, TableRow, Tableau},
    steps::{canonicalize, simplex_method},
};

// convert a system of equations into standard form
// returns (initial number of vars, number of slacks, number of artificials, variable_names, tableau)
// problem will have a primary objective with artificial variables
// and a secondary objective variable with the original objective
fn to_standard_form(mut problem: Problem) -> (usize, usize, usize, Vec<String>, Tableau) {
    // first add the slack and artificial variables
    let initial_vars = problem.num_vars();
    for constraint in &mut problem.constraints {
        // swap sign if necessary
        if constraint.row.rhs < Fraction::zero() {
            constraint.row *= -Fraction::one();
            constraint.inequality = constraint.inequality.opposite();
        }
    }
    // then count how many slacks and artificials we need;
    // do this before adding them that way we can put all artificials at the end
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

    // then add the slacks and artificials
    let mut added_slacks = 0;
    let mut added_artificials = 0;
    let mut tableau_rows = Vec::new();
    for constraint in &mut problem.constraints {
        // add extras
        for _ in 0..(total_slacks + total_artificial) {
            constraint.row.coefficients.push(Fraction::zero());
        }
        // then maybe add slack or artificial or both
        // and have whatever we use be the last variable in the row
        let basic_var = match constraint.inequality {
            Inequality::Le => {
                // since we have ... <= b, just do
                // ... + s = b
                constraint.row.coefficients[initial_vars + added_slacks] = Fraction::one();
                added_slacks += 1;
                initial_vars + added_slacks - 1
            }
            Inequality::Eq => {
                // since we have ... = b, do
                // ... + a = b, and we'd expect a to become 0
                constraint.row.coefficients[initial_vars + total_slacks + added_artificials] =
                    Fraction::one();
                added_artificials += 1;
                initial_vars + total_slacks + added_artificials - 1
            }
            Inequality::Ge => {
                // we have ... >= b, so we make it
                // ... - s = b, but then we also insert an artificial to use
                // as a pivot, giving
                // ... - s + a = b
                // slack first
                constraint.row.coefficients[initial_vars + added_slacks] = -Fraction::one();
                added_slacks += 1;
                // then artificial since that's what we'll use as our BV
                constraint.row.coefficients[initial_vars + total_slacks + added_artificials] =
                    Fraction::one();
                added_artificials += 1;
                initial_vars + total_slacks + added_artificials - 1
            }
        };
        tableau_rows.push(TableRow {
            row: constraint.row.clone(),
            basic_var,
        });
    }

    // then construct the objective row
    let mut objective_coeffs =
        vec![Fraction::zero(); initial_vars + total_slacks + total_artificial];
    for i in 0..total_artificial {
        objective_coeffs[initial_vars + total_slacks + i] = Fraction::one();
    }
    let mut original_objective = problem.objective.clone();
    original_objective
        .coefficients
        .extend(vec![Fraction::zero(); total_slacks + total_artificial]);
    let artificial_tableau = Tableau {
        objective: Row {
            coefficients: objective_coeffs,
            rhs: Fraction::zero(),
        },
        constraints: tableau_rows,
        secondary_objective: Some(original_objective),
    };

    let mut variable_names: Vec<_> = (0..initial_vars).map(|v| format!("x{}", v + 1)).collect();
    for i in 1..=total_slacks {
        variable_names.push(format!("s{i}"));
    }
    for i in 1..=total_artificial {
        variable_names.push(format!("a{i}"));
    }
    (
        initial_vars,
        total_slacks,
        total_artificial,
        variable_names,
        artificial_tableau,
    )
}

/// This function will attempt to solve a linear program
/// by running the two-pass simplex method:
/// The first pass finds a BFS using artificial variables, and the second pass finds a solution
/// using the BFS from the first pass.
///
/// If there is no solution or the objective is unbounded, it will error
/// All RHS values must be >= 0
///
/// Returns (first_pass_steps, second_pass_steps, variable_names)
pub fn solve(problem: Problem) -> Result<(Vec<Tableau>, Vec<Tableau>, Vec<String>), &'static str> {
    let (initial_vars, total_slacks, _, variable_names, artificial_tableau) =
        to_standard_form(problem);
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
            if !removed_row[row_idx] && tableau_row.basic_var >= initial_vars + total_slacks {
                Some(row_idx)
            } else {
                None
            }
        })
        .next()
    {
        let pivot_var = (0..(initial_vars + total_slacks))
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
                    .truncate(initial_vars + total_slacks);
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
        .truncate(initial_vars + total_slacks);

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

#[cfg(test)]
mod tests {
    use fraction::{Fraction, Zero};

    use super::solve;
    use crate::problem::{Inequality, InputConstraint, Problem, Row, Tableau};

    use Inequality::{Eq, Ge, Le};

    type Constraint<'a> = (&'a [i64], Inequality, i64);

    fn frac(x: i64) -> Fraction {
        Fraction::from(x)
    }

    fn row(coefficients: &[i64], rhs: i64) -> Row {
        Row {
            coefficients: coefficients.iter().map(|&c| frac(c)).collect(),
            rhs: frac(rhs),
        }
    }

    fn problem(constraints: &[Constraint], objective: &[i64]) -> Problem {
        Problem::new(
            constraints
                .iter()
                .map(|&(coefficients, inequality, rhs)| InputConstraint {
                    row: row(coefficients, rhs),
                    inequality,
                })
                .collect(),
            row(objective, 0),
        )
    }

    /// Values of every column in the tableau's basic solution
    fn basic_solution(tableau: &Tableau) -> Vec<Fraction> {
        let mut values = vec![Fraction::zero(); tableau.objective.coefficients.len()];
        for table_row in &tableau.constraints {
            values[table_row.basic_var] = table_row.row.rhs;
        }
        values
    }

    /// Solves the problem, checks that the answer is a feasible point with the expected
    /// objective value, and returns the values of the original variables
    fn assert_optimal(
        constraints: &[Constraint],
        objective: &[i64],
        expected: i64,
    ) -> Vec<Fraction> {
        let (first_pass, second_pass, names) =
            solve(problem(constraints, objective)).expect("expected an optimal solution");

        // every step of both phases should be a basic feasible solution
        for tableau in first_pass.iter().chain(&second_pass) {
            for table_row in &tableau.constraints {
                assert!(table_row.row.rhs >= Fraction::zero(), "RHS went negative");
            }
        }
        // phase 2 shouldn't have any artificial columns
        let last = second_pass.last().unwrap();
        for name in &names[..last.objective.coefficients.len()] {
            assert!(!name.starts_with('a'), "artificial {name} in phase 2");
        }

        let values = basic_solution(last);
        let x = &values[..objective.len()];
        for &(coefficients, inequality, rhs) in constraints {
            let lhs: Fraction = coefficients.iter().zip(x).map(|(&c, v)| frac(c) * *v).sum();
            let satisfied = match inequality {
                Le => lhs <= frac(rhs),
                Eq => lhs == frac(rhs),
                Ge => lhs >= frac(rhs),
            };
            assert!(satisfied, "solution violates {coefficients:?} {rhs}");
        }
        let value: Fraction = objective.iter().zip(x).map(|(&c, v)| frac(c) * *v).sum();
        assert_eq!(value, frac(expected), "objective value of the solution");
        assert_eq!(
            -last.objective.rhs,
            frac(expected),
            "objective value in the tableau"
        );
        x.to_vec()
    }

    fn assert_error(constraints: &[Constraint], objective: &[i64], expected: &str) {
        match solve(problem(constraints, objective)) {
            Ok(_) => panic!("expected error containing {expected:?}"),
            Err(e) => assert!(e.contains(expected), "got error {e:?}"),
        }
    }

    #[test]
    fn only_le_constraints() {
        // no artificials needed
        assert_optimal(&[(&[1, 1], Le, 4)], &[-1, -1], -4);
    }

    #[test]
    fn textbook_maximization() {
        // max 3x1 + 5x2 is min -3x1 - 5x2; optimum at (2, 6)
        let x = assert_optimal(
            &[(&[1, 0], Le, 4), (&[0, 2], Le, 12), (&[3, 2], Le, 18)],
            &[-3, -5],
            -36,
        );
        assert_eq!(x, vec![frac(2), frac(6)]);
    }

    #[test]
    fn ge_constraint() {
        assert_optimal(&[(&[1, 1], Ge, 2)], &[1, 1], 2);
    }

    #[test]
    fn eq_and_le_constraints() {
        let x = assert_optimal(&[(&[1, 1], Eq, 3), (&[1, 0], Le, 1)], &[1, 2], 5);
        assert_eq!(x, vec![frac(1), frac(2)]);
    }

    #[test]
    fn negative_rhs_le_becomes_ge() {
        // x1 - x2 <= -1 needs an artificial once it's flipped to -x1 + x2 >= 1
        assert_optimal(&[(&[1, -1], Le, -1)], &[0, 1], 1);
    }

    #[test]
    fn negative_rhs_ge_becomes_le() {
        // -x1 >= -3 is x1 <= 3
        assert_optimal(&[(&[-1], Ge, -3)], &[-1], -3);
    }

    #[test]
    fn redundant_equality_is_removed() {
        assert_optimal(&[(&[1, 1], Eq, 2), (&[1, 1], Eq, 2)], &[1, 3], 2);
    }

    #[test]
    fn artificial_left_at_zero_is_pivoted_out() {
        // phase 1 ends immediately with both artificials basic at 0;
        // the first is pivoted out and the second row is redundant
        assert_optimal(
            &[(&[1, -1], Eq, 0), (&[-1, 1], Eq, 0), (&[1, 0], Le, 3)],
            &[-1, -1],
            -6,
        );
    }

    #[test]
    fn all_rows_redundant() {
        assert_optimal(&[(&[0, 0], Eq, 0)], &[1, 0], 0);
    }

    #[test]
    fn all_rows_redundant_unbounded() {
        assert_error(&[(&[0, 0], Eq, 0)], &[-1, 0], "nbounded");
    }

    #[test]
    fn degenerate_beale_example() {
        // Beale's example, which cycles under some pivoting rules
        // (rows and objective scaled to integers)
        assert_optimal(
            &[
                (&[1, -32, -4, 36], Le, 0),
                (&[1, -24, -1, 6], Le, 0),
                (&[0, 0, 1, 0], Le, 1),
            ],
            &[-3, 80, -2, 24],
            -5,
        );
    }

    #[test]
    fn infeasible_inequalities() {
        assert_error(
            &[(&[1, 1], Le, 1), (&[1, 1], Ge, 2)],
            &[1, 1],
            "No feasible solution",
        );
    }

    #[test]
    fn infeasible_equality() {
        // x1 = -1 can't hold since x1 >= 0
        assert_error(&[(&[1], Eq, -1)], &[1], "No feasible solution");
    }

    #[test]
    fn unbounded() {
        assert_error(&[(&[1, -1], Ge, 1)], &[-1, 0], "nbounded");
    }

    #[test]
    fn variable_names() {
        let (_, _, names) = solve(problem(
            &[(&[1, 1], Le, 4), (&[1, 0], Ge, 1), (&[0, 1], Eq, 1)],
            &[1, 1],
        ))
        .unwrap();
        assert_eq!(names, ["x1", "x2", "s1", "s2", "a1", "a2"]);
    }
}
