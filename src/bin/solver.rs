extern crate Rubiks2x2x2_OptimalSolver;
use Rubiks2x2x2_OptimalSolver::solver;

fn main() {
    println!("starting");
    let solution = solver::solve("UUUURRRRFFFFDDDDLLLLBBBB");
    println!("{:?}", solution);
    println!("ok");
}
