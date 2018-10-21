extern crate Rubiks2x2x2_OptimalSolver;
use Rubiks2x2x2_OptimalSolver::solver;

fn main() {
    println!("starting");
    //let solution = solver::solve("UUUURRRRFFFFDDDDLLLLBBBB"); // normal cube
    //let solution = solver::solve("UUUURRFFFFLLDDDDLLBBBBRR"); // D
    let solution = solver::solve("RULBUFDBULRFBLDUDFLDRFRB"); // U R2 F' U' F R U R' F' U
    //let solution = solver::solve("BUUURRRRFFFFDDDDLLLLUBBB"); // invalid cube

    // U = White
    // F = Green
    // D = Yellow
    // B = Blue
    // L = Orange
    // R = Red
    println!("{:?}", solution);
    println!("ok");
}
