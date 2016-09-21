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

void printResidualError(std::string textBefore, const MatrixXd &A, const MatrixXd &xApproximate, const VectorXd &b)
{

    // This function computes the residual error ||b - Ax|| given the 3 inputs
    double residualError = (b - A*xApproximate).norm();
    // Print out the residual error in scientific notation
    std::cout << textBefore << std::endl
              << "Residual Error" << std::scientific << residualError << std::endl;
    //restore back the stream to fixed
    std::cout << std::fixed;

}

double maxApproximationError(const VectorXd &xApprox, const VectorXd &uExact)
{
    VectorXd diffVector = xApprox - uExact;
    // Get the absolute maximum coefficent of the error
    return diffVector.cwiseAbs().maxCoeff();
}


VectorXd residual_based_solver(const MatrixXd &A, const VectorXd &b, const VectorXd &x0, double tol, IterationMethod &iterMethod, double w) {
    //w is only used when it is referred to. Otherwise this value does not matter
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


    VectorXd b_new = iterMethod.calculateBnew(lower_A, upper_A, diagonal_A, D_inverse, b, x0, w); // calculate b_new

    int ic = 0; // iteration count

    while (residual.norm() > stopIterResidual) {
        x = iterMethod.calculateIteration(lower_A, upper_A, diagonal_A, D_inverse, b_new, x, w);
        residual = b - copiedA * x;
        ic += 1;
        std::string toPrint = "After iteration ";
        toPrint += std::to_string(ic);
        toPrint += ":";
        printCSVMatrix(toPrint, x);
    }
//    std::cout << ic;

    return x;
}

VectorXd consecutive_approximation_solver(const MatrixXd &A, const VectorXd &b, const VectorXd &x0, double tol, IterationMethod &iterMethod, double w) {
    //w is only used when it is referred to. Otherwise this value does not matter
    MatrixXd copiedA(A);
    VectorXd xOld(x0);
    VectorXd xNew;
    VectorXd diff(x0);
    diff.setOnes(); //initalize it to all ones

    // Input splitting
    MatrixXd lower_A = copiedA.triangularView<Eigen::StrictlyLower>();
    MatrixXd upper_A = copiedA.triangularView<Eigen::StrictlyUpper>();
    MatrixXd diagonal_A = copiedA.diagonal().asDiagonal();

    // Preparation for calculation
    MatrixXd D_inverse = diagonal_A.inverse();
    VectorXd b_new = iterMethod.calculateBnew(lower_A, upper_A, diagonal_A, D_inverse, b, x0, w); // calculate b_new

    int ic = 0; // iteration count

    while (diff.norm() > tol) {
        xNew = iterMethod.calculateIteration(lower_A, upper_A, diagonal_A, D_inverse, b_new, xOld, w);
        diff = xNew - xOld;
        xOld = xNew;

        ic += 1;
    }
    std::cout << "Iter count: " << ic << std::endl;

    return xNew;
}


#endif //CPPCODETEST_STOPPINGCRITERIONSOLVERS_HPP
