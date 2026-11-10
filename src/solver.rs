/// This function will attempt to solve a linear program
/// by running the simplex method twice:
/// once to find a BFS, once to find a solution.
/// If there is no solution or the objective is unbounded, it will panic
/// All RHS values must be >= 0
pub fn solve(equations: Vec<(Vec<Fraction>, Fraction)>, objective: Vec<Fraction>) -> Vec<Step> {
    // first, artificially augment the equations to force each one
    // to have an artificial variable
    // and then we run solve_with_bfs on that
    if equations.iter().any(|f| f.1 < Fraction::zero()) {
        panic!("RHS must be >= 0");
    }
    // now add artificial variables
    let original_num_vars = equations[0].0.len();
    let mut new_equations = equations.clone();
    let mut new_objective = vec![Fraction::zero(); original_num_vars];
    let mut new_vars = Vec::new();
    for i in 0..equations.len() {
        for j in 0..equations.len() {
            let val = if i == j {
                Fraction::one()
            } else {
                Fraction::zero()
            };
            new_equations[i].0.push(val);
        }
        // we want to minimize the artificial variables
        new_objective.push(Fraction::one());
        new_vars.push(i + original_num_vars);
    }
    // and run solve_with_bfs
    let first_part_steps = solve_with_bfs(new_equations, new_objective, new_vars);
    // we know that it's possible iff the objective was minimized to 0
    if first_part_steps.last().unwrap().objective.1 != Fraction::zero() {
        panic!("There are no solutions to this problem");
    }
    // now we just need to see what variables we used
    let mut used_rows = vec![false; equations.len()];
    let mut used_variable = vec![false; original_num_vars];
    let mut variables = Vec::new();
    let mut keep_row = vec![true; equations.len()];
    for step in first_part_steps.iter().rev() {
        if !used_rows[step.pivot_row] {
            if step.pivot_var < original_num_vars {
                assert!(keep_row[step.pivot_row]);
                used_rows[step.pivot_row] = true;
                used_variable[step.pivot_var] = true;
                variables.push(step.pivot_var);
            } else {
                // if we chose to use an artificial variable, it must be because there's
                // nothing else in the row
                assert!(
                    first_part_steps.last().unwrap().equations[step.pivot_row].1
                        == Fraction::zero()
                );
                // it doesn't matter what we choose; just take a variable with nonzero coefficient
                let mut nonzero_coef = None;
                for i in 0..original_num_vars {
                    if !used_variable[i]
                        && first_part_steps.last().unwrap().equations[step.pivot_row].0[i]
                            != Fraction::zero()
                    {
                        nonzero_coef = Some(i);
                        break;
                    }
                }
                if let Some(i) = nonzero_coef {
                    // otherwise, we can use it
                    used_rows[step.pivot_row] = true;
                    used_variable[i] = true;
                    variables.push(i);
                } else {
                    // drop the row since it's useless
                    keep_row[step.pivot_row] = false;
                }
            }
        }
    }
    assert!(variables.len() == keep_row.iter().filter(|x| **x).count());
    // then filter
    let mut new_equations = Vec::new();
    for (row, keep) in equations.iter().zip(keep_row.iter()) {
        if *keep {
            new_equations.push(row.clone());
        }
    }
    if new_equations.len() == 0 {
        if objective.iter().all(|x| *x > Fraction::zero()) {
            // one solution: have all variables = 0
            return vec![Step {
                equations: vec![(vec![Fraction::zero(); original_num_vars], Fraction::zero())],
                objective: (objective, Fraction::zero()),
                pivot_row: 0,
                pivot_var: 0,
            }];
        } else {
            // infinitely many solutions since there's at least one variable
            // with objective value <= 0. Then, we can just take any nonnegative real value for that
            // to get infinite solutions
            panic!("Infinitely many solutions");
        }
    }
    // now we just need to find a BFS
    // and run solve_with_bfs again
    solve_with_bfs(new_equations, objective, variables)
}
