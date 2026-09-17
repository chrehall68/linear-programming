mod lp;

use fraction::{Fraction, Zero};
use lp::{canonicalize, solve_with_bfs};
use std::io::{self, Write};
use std::panic;

use crate::lp::solve;

#[derive(Clone, Copy)]
enum Mode {
    Canonicalize,
    SolveWithBfs,
    Solve,
}

fn read_mode() -> Mode {
    println!("\nChoose a mode:");
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

/// Extracts a human-readable message from a caught panic payload, falling
/// back to a generic message when the payload isn't a plain string (as
/// produced by `panic!`/`assert!`).
fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        s.to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown error".to_string()
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

fn print_tableau(
    title: &str,
    equations: &[(Vec<Fraction>, Fraction)],
    objective: &(Vec<Fraction>, Fraction),
    num_vars: usize,
    basic_var_for_row: &[Option<usize>],
) {
    println!("\n{title}");

    let mut headers = vec!["basic".to_string()];
    headers.extend((0..num_vars).map(|v| format!("x{}", v + 1)));
    headers.push("RHS".to_string());

    let mut rows = Vec::new();
    for (i, (coeffs, rhs)) in equations.iter().enumerate() {
        let label = match basic_var_for_row[i] {
            Some(v) => format!("x{}", v + 1),
            None => format!("R{}", i + 1),
        };
        let mut row = vec![label];
        row.extend(coeffs.iter().map(format_fraction));
        row.push(format_fraction(rhs));
        rows.push(row);
    }

    let mut obj_row = vec!["obj".to_string()];
    obj_row.extend(objective.0.iter().map(format_fraction));
    obj_row.push(format_fraction(&objective.1));
    let obj_index = rows.len();
    rows.push(obj_row);

    print_table(&headers, &rows, &[obj_index]);
}

fn main() {
    println!("=== Linear Program Canonicalization ===");
    println!("Enter a system of equations in standard form (Ax = b) plus an objective row,");
    println!("then choose which variables should be basic in your basic feasible solution.\n");

    let num_vars = read_usize("Number of variables: ");
    let num_constraints = loop {
        let m = read_usize("Number of constraints: ");
        if m <= num_vars {
            break m;
        }
        println!("  → number of constraints can't exceed number of variables ({num_vars})");
    };

    println!(
        "\nFor each constraint, enter {num_vars} coefficient(s) then the RHS, space-separated."
    );
    println!("Example (n=3): 2 1 -1 10   means  2x1 + x2 - x3 = 10\n");

    let mut equations = Vec::with_capacity(num_constraints);
    for i in 0..num_constraints {
        let row = read_fraction_row(&format!("Constraint {}: ", i + 1), num_vars + 1);
        let rhs = row[num_vars];
        let coeffs = row[..num_vars].to_vec();
        equations.push((coeffs, rhs));
    }

    println!("\nEnter the {num_vars} objective coefficient(s), as in z = c1*x1 + c2*x2 + ...");
    let objective = read_fraction_row("Objective: ", num_vars);

    let mode = read_mode();
    match mode {
        Mode::Canonicalize => println!(
            "\nPick {num_constraints} variable(s) (by number, 1-{num_vars}) to be basic, in the exact pivot order to use."
        ),
        Mode::SolveWithBfs => println!(
            "\nPick {num_constraints} variable(s) (by number, 1-{num_vars}) for the starting basic feasible solution."
        ),
        Mode::Solve => (),
    }

    // Set a silent panic hook while we probe candidate bases, since an
    // incompatible choice of basic variables makes canonicalize()/solve() assert/panic.
    let default_hook = panic::take_hook();

    let steps = loop {
        match mode {
            Mode::Solve => break solve(equations.clone(), objective.clone()),
            _ => {
                panic::set_hook(Box::new(|_| {}));

                let variables = read_index_row(
                    "Basic variables, space-separated: ",
                    num_constraints,
                    num_vars,
                );

                let eqs = equations.clone();
                let obj = objective.clone();
                let attempt = panic::catch_unwind(panic::AssertUnwindSafe(move || match mode {
                    Mode::Canonicalize => canonicalize(eqs, obj, variables),
                    Mode::SolveWithBfs => solve_with_bfs(eqs, obj, variables),
                    _ => unreachable!(),
                }));

                match attempt {
                    Ok(steps) => break steps,
                    Err(e) => println!(
                        "  → that didn't work ({}). Try different variables.",
                        panic_message(&*e)
                    ),
                }
            }
        }
    };
    panic::set_hook(default_hook);

    println!();
    let mut basic_var_for_row: Vec<Option<usize>> = vec![None; num_constraints];
    print_tableau(
        "Initial tableau",
        &equations,
        &(objective.clone(), Fraction::zero()),
        num_vars,
        &basic_var_for_row,
    );

    for (step_num, step) in steps.iter().enumerate() {
        basic_var_for_row[step.pivot_row] = Some(step.pivot_var);

        print_tableau(
            &format!(
                "Step {}: pivot on x{} (row {})",
                step_num + 1,
                step.pivot_var + 1,
                step.pivot_row + 1
            ),
            &step.equations,
            &step.objective,
            num_vars,
            &basic_var_for_row,
        );
    }

    println!("\n=== Basic Feasible Solution ===");
    let final_step = steps.last().unwrap();
    let mut values = vec![Fraction::zero(); num_vars];
    for (row, var) in basic_var_for_row.iter().enumerate() {
        if let Some(v) = *var {
            values[v] = final_step.equations[row].1;
        }
    }
    for (i, v) in values.iter().enumerate() {
        println!("  x{} = {}", i + 1, v);
    }
    // The objective row started life as "c·x = 0"; after fully pivoting it (like the
    // other rows) its reduced coefficients are zero on every basic variable, so at the
    // BFS (non-basic vars = 0) the row reads "0 = <final RHS>" only because the true
    // objective value is the *negative* of that final RHS.
    println!("  objective value = {}", -final_step.objective.1);
}
