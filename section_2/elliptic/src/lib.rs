//! This crate provides a means of comparing the convergence of several different schemes for the Laplace's equations.
//!
//! Section 2.4 of the book introduces the relaxation methods for the Laplace's equations and discusses the convergence of them.
//!
//! All of the methods mentioned in the book are implemented in this crate.
//!
//! Using this crate, you can actually compute and see the convergence of each method.

pub mod input;
pub mod output;
pub mod solver;

use solver::Solver;
use std::error::Error;
use std::io::Write;

/// Run the solver and output the results.
pub fn run(solver: &mut impl Solver, outputstream: &mut impl Write) -> Result<(), Box<dyn Error>> {
    // calculate and output
    solver.exec()?;
    output::output(outputstream, solver.borrow_u())?;
    println!(
        "The solution is converged at {} iterations.",
        solver.get_n_iter()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::prelude::*;
    use solver::point_jacobi_solver::{PointJacobiSolver, PointJacobiSolverNewParams};
    use solver::sor_solver::{SorSolver, SorSolverNewParams};

    #[test]
    fn fn_run_works_with_point_jacobi_solver() {
        // setup input parameters
        let n_x = 8;
        let n_y = 8;

        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup initial and boundary conditions
        let mut u_init: Array2<f64> = Array::zeros((n_x + 1, n_y + 1));
        u_init.slice_mut(s![.., n_y]).assign(&Array::ones(n_x + 1));

        // initialize the solver
        let new_params = PointJacobiSolverNewParams {
            u_init,
            n_iter_max: 300,
        };
        let mut solver = PointJacobiSolver::new(new_params).unwrap();

        // execute run()
        run(&mut solver, &mut outputstream).unwrap();

        // check if the output is correct
        let output_expected = "\
0 0 0.0000000000
0 1 0.0000000000
0 2 0.0000000000
0 3 0.0000000000
0 4 0.0000000000
0 5 0.0000000000
0 6 0.0000000000
0 7 0.0000000000
0 8 1.0000000000

1 0 0.0000000000
1 1 0.0174130457
1 2 0.0376833281
1 3 0.0643733244
1 4 0.1029411760
1 5 0.1635678513
1 6 0.2693019654
1 7 0.4825869540
1 8 1.0000000000

2 0 0.0000000000
2 1 0.0319688546
2 2 0.0689469426
2 3 0.1168687935
2 4 0.1838235286
2 5 0.2820282638
2 6 0.4310530562
2 7 0.6610458507
2 8 1.0000000000

3 0 0.0000000000
3 1 0.0415154301
3 2 0.0892667944
3 3 0.1503313786
3 4 0.2334558813
3 5 0.3496686194
3 6 0.5118361453
3 7 0.7305433926
3 8 1.0000000000

4 0 0.0000000000
4 1 0.0448260716
4 2 0.0962734264
4 3 0.1617340455
4 4 0.2499999988
4 5 0.3713541876
4 6 0.5360795131
4 7 0.7492915745
4 8 1.0000000000

5 0 0.0000000000
5 1 0.0415154301
5 2 0.0892667944
5 3 0.1503313786
5 4 0.2334558813
5 5 0.3496686194
5 6 0.5118361453
5 7 0.7305433926
5 8 1.0000000000

6 0 0.0000000000
6 1 0.0319688546
6 2 0.0689469426
6 3 0.1168687935
6 4 0.1838235286
6 5 0.2820282638
6 6 0.4310530562
6 7 0.6610458507
6 8 1.0000000000

7 0 0.0000000000
7 1 0.0174130457
7 2 0.0376833281
7 3 0.0643733244
7 4 0.1029411760
7 5 0.1635678513
7 6 0.2693019654
7 7 0.4825869540
7 8 1.0000000000

8 0 0.0000000000
8 1 0.0000000000
8 2 0.0000000000
8 3 0.0000000000
8 4 0.0000000000
8 5 0.0000000000
8 6 0.0000000000
8 7 0.0000000000
8 8 1.0000000000

";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }

    #[test]
    fn fn_run_works_with_sor_solver() {
        // setup input parameters
        let n_x = 8;
        let n_y = 8;

        // setup output stream
        let mut outputstream: Vec<u8> = Vec::new();

        // setup initial and boundary conditions
        let mut u_init: Array2<f64> = Array::zeros((n_x + 1, n_y + 1));
        u_init.slice_mut(s![.., n_y]).assign(&Array::ones(n_x + 1));

        // initialize the solver
        let new_params = SorSolverNewParams {
            u_init,
            n_iter_max: 300,
            omega: 1.5,
        };
        let mut solver = SorSolver::new(new_params).unwrap();

        // execute run()
        run(&mut solver, &mut outputstream).unwrap();

        // check if the output is correct
        let output_expected = "\
0 0 0.0000000000
0 1 0.0000000000
0 2 0.0000000000
0 3 0.0000000000
0 4 0.0000000000
0 5 0.0000000000
0 6 0.0000000000
0 7 0.0000000000
0 8 1.0000000000

1 0 0.0000000000
1 1 0.0174130458
1 2 0.0376833284
1 3 0.0643733247
1 4 0.1029411764
1 5 0.1635678517
1 6 0.2693019657
1 7 0.4825869542
1 8 1.0000000000

2 0 0.0000000000
2 1 0.0319688548
2 2 0.0689469431
2 3 0.1168687942
2 4 0.1838235294
2 5 0.2820282646
2 6 0.4310530568
2 7 0.6610458510
2 8 1.0000000000

3 0 0.0000000000
3 1 0.0415154305
3 2 0.0892667951
3 3 0.1503313795
3 4 0.2334558823
3 5 0.3496686204
3 6 0.5118361460
3 7 0.7305433930
3 8 1.0000000000

4 0 0.0000000000
4 1 0.0448260720
4 2 0.0962734272
4 3 0.1617340466
4 4 0.2500000000
4 5 0.3713541887
4 6 0.5360795139
4 7 0.7492915750
4 8 1.0000000000

5 0 0.0000000000
5 1 0.0415154305
5 2 0.0892667951
5 3 0.1503313796
5 4 0.2334558823
5 5 0.3496686204
5 6 0.5118361460
5 7 0.7305433930
5 8 1.0000000000

6 0 0.0000000000
6 1 0.0319688549
6 2 0.0689469432
6 3 0.1168687942
6 4 0.1838235294
6 5 0.2820282646
6 6 0.4310530568
6 7 0.6610458510
6 8 1.0000000000

7 0 0.0000000000
7 1 0.0174130458
7 2 0.0376833284
7 3 0.0643733248
7 4 0.1029411765
7 5 0.1635678517
7 6 0.2693019657
7 7 0.4825869542
7 8 1.0000000000

8 0 0.0000000000
8 1 0.0000000000
8 2 0.0000000000
8 3 0.0000000000
8 4 0.0000000000
8 5 0.0000000000
8 6 0.0000000000
8 7 0.0000000000
8 8 1.0000000000

";
        assert_eq!(String::from_utf8(outputstream).unwrap(), output_expected);
    }
}
