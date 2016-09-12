//
// Created by avi on 9/11/16.
//

#ifndef CPPCODETEST_ITERATIVEMETHODS_HPP
#define CPPCODETEST_ITERATIVEMETHODS_HPP
#include "Eigen/Dense"
#include <tuple>
#include <iostream>
#include <vector>

using namespace Eigen;

VectorXd Jacobi_Iteration(const MatrixXd &A, const VectorXd &b, VectorXd x0, double tol)
{
    MatrixXd copiedA(A);
    VectorXd x = x0;
    VectorXd r0 = b - copiedA*x;

    VectorXd residual = r0; // start at current residual
    double stopIterResidual = tol*residual.norm();

    MatrixXd lower_A = copiedA.triangularView<Eigen::StrictlyLower>();
    MatrixXd upper_A = copiedA.triangularView<Eigen::StrictlyUpper>();
    MatrixXd diagonal_A = copiedA.diagonal().asDiagonal();
//    std::cout << "Original:" <<std::endl << copiedA << std::endl;
//    std::cout << "Lower:" << std::endl << lower_A  << std::endl;
//    std::cout << "Upper" << std::endl << upper_A << std::endl;
//    std::cout << "Diagonal" << std::endl << diagonal_A << std::endl;
//    std::cout << lower_A + upper_A + diagonal_A - copiedA << std::endl;

    return x;

}




#endif //CPPCODETEST_ITERATIVEMETHODS_HPP
