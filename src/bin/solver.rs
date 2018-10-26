extern crate Rubiks2x2x2_OptimalSolver;
use Rubiks2x2x2_OptimalSolver as solver;
use std::env;

fn main() {
    println!("starting");

    let args: Vec<String> = env::args().collect();
    let solution;
    if args.len() < 2 {
        //solution = solver::solve("UUUURRRRFFFFDDDDLLLLBBBB"); // normal cube
        //solution = solver::solve("UUUURRFFFFLLDDDDLLBBBBRR"); // D
        solution = solver::solve("RULBUFDBULRFBLDUDFLDRFRB"); // U R2 F' U' F R U R' F' U
        //solution = solver::solve("BUUURRRRFFFFDDDDLLLLUBBB"); // invalid cube
    } else {
        let cubestring = &args[1];
        solution = solver::solve(cubestring);
    }

    // U = White
    // F = Green
    // D = Yellow
    // B = Blue
    // L = Orange
    // R = Red
    match solution {
        Ok(s) => { print!("{}", s); println!("ok"); },
        Err(s) => eprintln!("{}", s),
    }
}
