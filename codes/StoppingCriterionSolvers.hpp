//
// Created by avi on 9/12/16.
//

#ifndef CPPCODETEST_STOPPINGCRITERIONSOLVERS_HPP
#define CPPCODETEST_STOPPINGCRITERIONSOLVERS_HPP

#include "Eigen/Dense"
#include <tuple>
#include <iostream>
#include <vector>

#include "Solvers.hpp"
#include "IterativeMethods.hpp"

using namespace Eigen;


VectorXd residual_based_solver(const MatrixXd &A, const VectorXd &b, const VectorXd &x0, double tol, IterationMethod &iterMethod)
{
    MatrixXd copiedA(A);
    VectorXd x = x0;
    VectorXd r0 = b - copiedA*x;

    VectorXd residual = r0; // start at current residual
    double stopIterResidual = tol*residual.norm();

    // Input splitting
    MatrixXd lower_A = copiedA.triangularView<Eigen::StrictlyLower>();
    MatrixXd upper_A = copiedA.triangularView<Eigen::StrictlyUpper>();
    MatrixXd diagonal_A = copiedA.diagonal().asDiagonal();

    // Preparation for calculation
    MatrixXd D_inverse = diagonal_A.inverse();

//    virtual VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd D_inverse, VectorXd b, VectorXd x)=0;

    VectorXd b_new = iterMethod.calculateBnew(lower_A, upper_A, D_inverse, b, x0 ); // calculate b_new

    int ic = 0; // iteration count

    while (residual.norm() > stopIterResidual) {
        x = iterMethod.calculateIteration(lower_A, upper_A, D_inverse, b_new, x);
        residual = b - copiedA * x;
        ic += 1;
    }
    std::cout << "Iter count: " << ic << std::endl;

    return x;
}


#endif //CPPCODETEST_STOPPINGCRITERIONSOLVERS_HPP
