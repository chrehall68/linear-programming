// convert a system of equations into canonical form
// for a given number of variables
// returns the intermediate steps too
use crate::problem::{Problem, Tableau};
use fraction::{Fraction, One, Zero};

/// Convert a system of equations into canonical form
/// requires that the assigned basic variables are linearly independent
/// Returns [initial tableau, intermediate steps]. If the tableau is already in canonical form,
/// only returns [initial tableau]
fn canonicalize(mut prev_tableau: Tableau) -> Vec<Tableau> {
    let mut steps = vec![prev_tableau.clone()];
    for pivot_row in 0..prev_tableau.constraints.len() {
        // find a row where the chosen variable is not zero
        // and we haven't already used that row as a pivot
        let tableau_row = &prev_tableau.constraints[pivot_row];
        let pivot_var = tableau_row.basic_var;
        let pivot_val = tableau_row.row.coefficients[pivot_var];
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
    steps
}

/// attempt to solve a linear program
/// given an equation in standard form and a set of variables
/// to start with.
/// The given variables must form a basic feasible solution,
/// and if they don't, this will panic
pub fn solve_with_bfs(initial_tableau: Tableau) -> Vec<Tableau> {
    let num_vars = initial_tableau.constraints[0].row.coefficients.len();
    let mut steps = canonicalize(initial_tableau);
    let mut last = steps.last().unwrap().clone();

    // first make sure we're in a feasible solution
    // if we are, then all other things we visit will also be feasible
    for table_row in &last.constraints {
        if table_row.row.rhs < Fraction::zero() {
            panic!("not in a feasible solution");
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
        for var in 0..num_vars {
            if last.objective.coefficients[var] < Fraction::zero() {
                // make sure that the variable is bounded
                let bounded = last
                    .constraints
                    .iter()
                    .any(|table_row| table_row.row.coefficients[var] > Fraction::zero());
                if !bounded {
                    panic!("objective is unbounded");
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
        let mut next_steps = canonicalize(last.clone());
        steps.append(&mut next_steps);
        last = steps.last().unwrap().clone();
    }
    steps
}

// /// This function will attempt to solve a linear program
// /// by running the simplex method twice:
// /// once to find a BFS, once to find a solution.
// /// If there is no solution or the objective is unbounded, it will panic
// /// All RHS values must be >= 0
// pub fn solve(equations: Vec<(Vec<Fraction>, Fraction)>, objective: Vec<Fraction>) -> Vec<Step> {
//     // first, artificially augment the equations to force each one
//     // to have an artificial variable
//     // and then we run solve_with_bfs on that
//     if equations.iter().any(|f| f.1 < Fraction::zero()) {
//         panic!("RHS must be >= 0");
//     }
//     // now add artificial variables
//     let original_num_vars = equations[0].0.len();
//     let mut new_equations = equations.clone();
//     let mut new_objective = vec![Fraction::zero(); original_num_vars];
//     let mut new_vars = Vec::new();
//     for i in 0..equations.len() {
//         for j in 0..equations.len() {
//             let val = if i == j {
//                 Fraction::one()
//             } else {
//                 Fraction::zero()
//             };
//             new_equations[i].0.push(val);
//         }
//         // we want to minimize the artificial variables
//         new_objective.push(Fraction::one());
//         new_vars.push(i + original_num_vars);
//     }
//     // and run solve_with_bfs
//     let first_part_steps = solve_with_bfs(new_equations, new_objective, new_vars);
//     // we know that it's possible iff the objective was minimized to 0
//     if first_part_steps.last().unwrap().objective.1 != Fraction::zero() {
//         panic!("There are no solutions to this problem");
//     }
//     // now we just need to see what variables we used
//     let mut used_rows = vec![false; equations.len()];
//     let mut used_variable = vec![false; original_num_vars];
//     let mut variables = Vec::new();
//     let mut keep_row = vec![true; equations.len()];
//     for step in first_part_steps.iter().rev() {
//         if !used_rows[step.pivot_row] {
//             if step.pivot_var < original_num_vars {
//                 assert!(keep_row[step.pivot_row]);
//                 used_rows[step.pivot_row] = true;
//                 used_variable[step.pivot_var] = true;
//                 variables.push(step.pivot_var);
//             } else {
//                 // if we chose to use an artificial variable, it must be because there's
//                 // nothing else in the row
//                 assert!(
//                     first_part_steps.last().unwrap().equations[step.pivot_row].1
//                         == Fraction::zero()
//                 );
//                 // it doesn't matter what we choose; just take a variable with nonzero coefficient
//                 let mut nonzero_coef = None;
//                 for i in 0..original_num_vars {
//                     if !used_variable[i]
//                         && first_part_steps.last().unwrap().equations[step.pivot_row].0[i]
//                             != Fraction::zero()
//                     {
//                         nonzero_coef = Some(i);
//                         break;
//                     }
//                 }
//                 if let Some(i) = nonzero_coef {
//                     // otherwise, we can use it
//                     used_rows[step.pivot_row] = true;
//                     used_variable[i] = true;
//                     variables.push(i);
//                 } else {
//                     // drop the row since it's useless
//                     keep_row[step.pivot_row] = false;
//                 }
//             }
//         }
//     }
//     assert!(variables.len() == keep_row.iter().filter(|x| **x).count());
//     // then filter
//     let mut new_equations = Vec::new();
//     for (row, keep) in equations.iter().zip(keep_row.iter()) {
//         if *keep {
//             new_equations.push(row.clone());
//         }
//     }
//     if new_equations.len() == 0 {
//         if objective.iter().all(|x| *x > Fraction::zero()) {
//             // one solution: have all variables = 0
//             return vec![Step {
//                 equations: vec![(vec![Fraction::zero(); original_num_vars], Fraction::zero())],
//                 objective: (objective, Fraction::zero()),
//                 pivot_row: 0,
//                 pivot_var: 0,
//             }];
//         } else {
//             // infinitely many solutions since there's at least one variable
//             // with objective value <= 0. Then, we can just take any nonnegative real value for that
//             // to get infinite solutions
//             panic!("Infinitely many solutions");
//         }
//     }
//     // now we just need to find a BFS
//     // and run solve_with_bfs again
//     solve_with_bfs(new_equations, objective, variables)
// }
