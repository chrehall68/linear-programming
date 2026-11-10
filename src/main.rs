mod problem;
mod solver;
mod steps;

use fraction::{Fraction, One, Zero};
use problem::{Inequality, InputConstraint, Problem, Row, TableRow, Tableau};
use std::io::{self, Write};
use steps::{canonicalize, simplex_method};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Canonicalize,
    SolveWithBfs,
    Solve,
}

fn read_mode() -> Mode {
    println!("Choose a mode:");
    println!("  1) canonicalize - reduce to canonical form using an exact pivot order you choose");
    println!(
        "  2) solve with bfs  - run the simplex method starting from a basic feasible solution you choose"
    );
    println!("  3) solve         - run the simplex method starting with artificial variables");
    loop {
        match prompt("Mode (1/2/3): ").as_str() {
            "1" => return Mode::Canonicalize,
            "2" => return Mode::SolveWithBfs,
            "3" => return Mode::Solve,
            _ => println!("  → please enter 1, 2, or 3"),
        }
    }
}

fn prompt(msg: &str) -> String {
    print!("{msg}");
    io::stdout().flush().unwrap();
    let mut line = String::new();
    let bytes_read = io::stdin()
        .read_line(&mut line)
        .expect("failed to read input");
    if bytes_read == 0 {
        // stdin closed (EOF) with no more input to give us
        println!("\n(input ended, exiting)");
        std::process::exit(0);
    }
    line.trim().to_string()
}

fn read_usize(msg: &str) -> usize {
    loop {
        match prompt(msg).parse::<usize>() {
            Ok(v) if v > 0 => return v,
            _ => println!("  → please enter a positive whole number"),
        }
    }
}

/// Reads a line of `count` whitespace-separated fractions (integers, decimals
/// like `0.5`, or fractions like `1/2`, negatives allowed).
fn read_fraction_row(msg: &str, count: usize) -> Vec<Fraction> {
    loop {
        let line = prompt(msg);
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() != count {
            println!("  → expected {count} value(s), got {}", tokens.len());
            continue;
        }
        let mut values = Vec::with_capacity(count);
        let mut ok = true;
        for t in &tokens {
            match t.parse::<Fraction>() {
                Ok(f) => values.push(f),
                Err(_) => {
                    println!("  → couldn't parse '{t}' as a number (try 3, -2, 1/2, 0.75)");
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return values;
        }
    }
}

/// Reads `count` distinct 1-indexed variable numbers in [1, max], returns them 0-indexed.
fn read_index_row(msg: &str, count: usize, max: usize) -> Vec<usize> {
    loop {
        let line = prompt(msg);
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() != count {
            println!("  → expected {count} index/indices, got {}", tokens.len());
            continue;
        }
        let mut values = Vec::with_capacity(count);
        let mut ok = true;
        for t in &tokens {
            match t.parse::<usize>() {
                Ok(v) if v >= 1 && v <= max => values.push(v - 1),
                _ => {
                    println!("  → '{t}' must be a whole number between 1 and {max}");
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue;
        }
        let mut sorted = values.clone();
        sorted.sort();
        sorted.dedup();
        if sorted.len() != values.len() {
            println!("  → indices must be distinct");
            continue;
        }
        return values;
    }
}

fn format_fraction(f: &Fraction) -> String {
    if f.is_zero() {
        // the fraction crate can print a negative-zero result as "-0"; normalize it
        "0".to_string()
    } else {
        format!("{f}")
    }
}

fn center(s: &str, width: usize) -> String {
    let len = s.chars().count();
    if len >= width {
        return s.to_string();
    }
    let total = width - len;
    let left = total / 2;
    let right = total - left;
    format!("{}{}{}", " ".repeat(left), s, " ".repeat(right))
}

/// Prints a box-drawn table. `separators_before` lists row indices (into
/// `rows`) that should get a divider printed above them (used to set the
/// objective row apart from the equation rows).
fn print_table(headers: &[String], rows: &[Vec<String>], separators_before: &[usize]) {
    let mut widths: Vec<usize> = headers.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (c, cell) in row.iter().enumerate() {
            widths[c] = widths[c].max(cell.chars().count());
        }
    }
    let seg: Vec<String> = widths.iter().map(|w| "─".repeat(w + 2)).collect();
    let top = format!("┌{}┐", seg.join("┬"));
    let mid = format!("├{}┤", seg.join("┼"));
    let bot = format!("└{}┘", seg.join("┴"));

    let fmt_row = |cells: &[String]| -> String {
        let parts: Vec<String> = cells
            .iter()
            .enumerate()
            .map(|(c, v)| format!(" {} ", center(v, widths[c])))
            .collect();
        format!("│{}│", parts.join("│"))
    };

    println!("{top}");
    println!("{}", fmt_row(headers));
    println!("{mid}");
    for (i, row) in rows.iter().enumerate() {
        if separators_before.contains(&i) {
            println!("{mid}");
        }
        println!("{}", fmt_row(row));
    }
    println!("{bot}");
}

/// Reads a constraint of the form `c1 c2 ... cn <op> b`, where `<op>` is one of
/// `<=`, `>=`, or `=`.
fn read_constraint(msg: &str, count: usize) -> InputConstraint {
    loop {
        let line = prompt(msg);
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() != count + 2 {
            println!(
                "  → expected {count} coefficient(s), a relation, and a RHS, got {} value(s)",
                tokens.len()
            );
            continue;
        }
        let inequality = match tokens[count] {
            "<=" => Inequality::Le,
            ">=" => Inequality::Ge,
            "=" => Inequality::Eq,
            t => {
                println!("  → '{t}' must be one of <=, >=, =");
                continue;
            }
        };
        let mut values = Vec::with_capacity(count + 1);
        for t in tokens[..count].iter().chain(&tokens[count + 1..]) {
            match t.parse::<Fraction>() {
                Ok(f) => values.push(f),
                Err(_) => {
                    println!("  → couldn't parse '{t}' as a number (try 3, -2, 1/2, 0.75)");
                    break;
                }
            }
        }
        if values.len() != count + 1 {
            continue;
        }
        let rhs = values.pop().unwrap();
        return InputConstraint {
            row: Row {
                coefficients: values,
                rhs,
            },
            inequality,
        };
    }
}

fn read_maximize() -> bool {
    loop {
        match prompt("Minimize or maximize? (min/max): ").as_str() {
            "min" => return false,
            "max" => return true,
            _ => println!("  → please enter min or max"),
        }
    }
}

/// Whether `var` is in canonical form as the basic variable of `row`:
/// its column is the unit vector for `row` and its objective coefficient is zero
fn is_canonical(tableau: &Tableau, var: usize, row: usize) -> bool {
    let unit_column = tableau.constraints.iter().enumerate().all(|(i, table_row)| {
        let expected = if i == row {
            Fraction::one()
        } else {
            Fraction::zero()
        };
        table_row.row.coefficients[var] == expected
    });
    unit_column && tableau.objective.coefficients[var] == Fraction::zero()
}

fn print_tableau(title: &str, tableau: &Tableau, names: &[String], objective_names: (&str, &str)) {
    println!("\n{title}");

    let num_cols = tableau.objective.coefficients.len();
    let mut headers = vec!["basic".to_string()];
    headers.extend(names[..num_cols].iter().cloned());
    headers.push("RHS".to_string());

    let mut rows = Vec::new();
    for table_row in &tableau.constraints {
        let mut row = vec![names[table_row.basic_var].clone()];
        row.extend(table_row.row.coefficients.iter().map(format_fraction));
        row.push(format_fraction(&table_row.row.rhs));
        rows.push(row);
    }
    let obj_index = rows.len();
    let objectives = [Some(&tableau.objective), tableau.secondary_objective.as_ref()];
    for (label, objective) in [objective_names.0, objective_names.1].iter().zip(objectives) {
        if let Some(objective) = objective {
            let mut row = vec![label.to_string()];
            row.extend(objective.coefficients.iter().map(format_fraction));
            row.push(format_fraction(&objective.rhs));
            rows.push(row);
        }
    }

    print_table(&headers, &rows, &[obj_index]);
}

/// Prints each step, skipping the steps where a new basic variable was
/// assigned but not yet pivoted on (those look identical to the step before them).
fn print_steps(phase: &str, steps: &[Tableau], names: &[String], objective_names: (&str, &str)) {
    print_tableau(
        &format!("{phase}: initial tableau"),
        &steps[0],
        names,
        objective_names,
    );
    let mut step_num = 0;
    for window in steps.windows(2) {
        let (prev, next) = (&window[0], &window[1]);
        // the pivot row is the one whose basic variable just became canonical
        let pivot = (0..next.constraints.len()).find(|&row| {
            let var = next.constraints[row].basic_var;
            !is_canonical(prev, var, row) && is_canonical(next, var, row)
        });
        if let Some(row) = pivot {
            step_num += 1;
            let var = next.constraints[row].basic_var;
            print_tableau(
                &format!(
                    "{phase}, step {step_num}: pivot on {} (row {})",
                    names[var],
                    row + 1
                ),
                next,
                names,
                objective_names,
            );
        }
    }
}

fn print_solution(tableau: &Tableau, names: &[String], objective_value: Fraction) {
    let num_cols = tableau.objective.coefficients.len();
    let mut values = vec![Fraction::zero(); num_cols];
    for table_row in &tableau.constraints {
        values[table_row.basic_var] = table_row.row.rhs;
    }
    for (name, v) in names.iter().zip(&values) {
        println!("  {name} = {}", format_fraction(v));
    }
    println!("  objective value = {}", format_fraction(&objective_value));
}

fn main() {
    println!("=== Linear Program Solver ===\n");
    let mode = read_mode();

    let num_vars = read_usize("\nNumber of variables: ");
    let num_constraints = loop {
        let m = read_usize("Number of constraints: ");
        // modes 1 and 2 need a distinct basic variable for each constraint
        if mode == Mode::Solve || m <= num_vars {
            break m;
        }
        println!("  → number of constraints can't exceed number of variables ({num_vars})");
    };
    let names: Vec<String> = (0..num_vars).map(|v| format!("x{}", v + 1)).collect();

    if mode == Mode::Solve {
        println!(
            "\nFor each constraint, enter {num_vars} coefficient(s), then <=, >=, or =, then the RHS, space-separated."
        );
        println!("Example (n=3): 2 1 -1 <= 10   means  2x1 + x2 - x3 <= 10\n");
        let constraints: Vec<_> = (0..num_constraints)
            .map(|i| read_constraint(&format!("Constraint {}: ", i + 1), num_vars))
            .collect();

        println!("\nEnter the {num_vars} objective coefficient(s), as in z = c1*x1 + c2*x2 + ...");
        let mut objective = Row {
            coefficients: read_fraction_row("Objective: ", num_vars),
            rhs: Fraction::zero(),
        };
        // the solver minimizes, so maximize z by minimizing -z
        let maximize = read_maximize();
        if maximize {
            objective *= -Fraction::one();
        }

        match solver::solve(Problem::new(constraints, objective)) {
            Ok((first_pass, second_pass, names)) => {
                print_steps("Phase 1", &first_pass, &names, ("w", "z"));
                print_steps("Phase 2", &second_pass, &names, ("z", ""));

                println!("\n=== Optimal Solution ===");
                // the objective row's RHS is -z at the end
                let last = second_pass.last().unwrap();
                let value = if maximize {
                    last.objective.rhs
                } else {
                    -last.objective.rhs
                };
                print_solution(last, &names, value);
            }
            Err(e) => println!("\nNo optimal solution: {e}"),
        }
        return;
    }

    println!(
        "\nFor each constraint, enter {num_vars} coefficient(s) then the RHS, space-separated."
    );
    println!("Example (n=3): 2 1 -1 10   means  2x1 + x2 - x3 = 10\n");
    let equations: Vec<Row> = (0..num_constraints)
        .map(|i| {
            let mut coefficients =
                read_fraction_row(&format!("Constraint {}: ", i + 1), num_vars + 1);
            let rhs = coefficients.pop().unwrap();
            Row { coefficients, rhs }
        })
        .collect();

    println!("\nEnter the {num_vars} objective coefficient(s), as in z = c1*x1 + c2*x2 + ...");
    let objective = Row {
        coefficients: read_fraction_row("Objective: ", num_vars),
        rhs: Fraction::zero(),
    };

    match mode {
        Mode::Canonicalize => println!(
            "\nPick {num_constraints} variable(s) (by number, 1-{num_vars}) to be basic, in the exact pivot order to use."
        ),
        _ => println!(
            "\nPick {num_constraints} variable(s) (by number, 1-{num_vars}) for the starting basic feasible solution."
        ),
    }

    let steps = loop {
        let variables = read_index_row(
            "Basic variables, space-separated: ",
            num_constraints,
            num_vars,
        );
        let tableau = Tableau {
            constraints: equations
                .iter()
                .zip(variables)
                .map(|(row, basic_var)| TableRow {
                    row: row.clone(),
                    basic_var,
                })
                .collect(),
            objective: objective.clone(),
            secondary_objective: None,
        };
        let attempt = match mode {
            Mode::Canonicalize => canonicalize(tableau),
            _ => simplex_method(tableau),
        };
        match attempt {
            Ok(steps) => break steps,
            Err(e) => println!("  → that didn't work ({e}). Try different variables."),
        }
    };

    let phase = match mode {
        Mode::Canonicalize => "Canonicalize",
        _ => "Simplex",
    };
    print_steps(phase, &steps, &names, ("z", ""));

    println!("\n=== Basic Feasible Solution ===");
    // The objective row started life as "c·x = 0"; after fully pivoting it (like the
    // other rows) its reduced coefficients are zero on every basic variable, so at the
    // BFS (non-basic vars = 0) the row reads "0 = <final RHS>" only because the true
    // objective value is the *negative* of that final RHS.
    let last = steps.last().unwrap();
    print_solution(last, &names, -last.objective.rhs);
}
