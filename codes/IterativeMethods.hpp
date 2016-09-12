//
// Created by avi on 9/11/16.
//

#ifndef CPPCODETEST_ITERATIVEMETHODS_HPP
#define CPPCODETEST_ITERATIVEMETHODS_HPP
#include "Eigen/Dense"
#include <tuple>
#include <iostream>
#include <vector>

#include "Solvers.hpp"

using namespace Eigen;



class IterationMethod{
    public:
        virtual VectorXd
        calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd Dinverse, VectorXd b, VectorXd x0,
                              double w)=0;
        virtual VectorXd
        calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd Dinverse, VectorXd bNew, VectorXd x,
                                   double w)=0;

};

class JacobiIteration: public IterationMethod
{
public:
    VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd D_inverse, VectorXd b, VectorXd x0) {
        return D_inverse*b;
    }

    VectorXd calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd DInverse, VectorXd bNew, VectorXd x){
        return -1* DInverse * (lowerA * x + upperA * x) + bNew;
    }
};

class GaussSiedelIteration: public IterationMethod
{
public:
    VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd Dinverse, VectorXd b, VectorXd x0,
                   double w){

        return forward_subst(diagonalA + lowerA, b); // calculate b_new;
    }

    VectorXd calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd DInverse, VectorXd bNew, VectorXd x, double w) {
        return -1*forward_subst(diagonalA + lowerA, upperA*x) + bNew;
    }
};


class SORIteration: public IterationMethod
{
public:
    VectorXd
    calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd D_inverse, VectorXd b, VectorXd x0, double w) {
        return w*forward_subst(diagonalA + w*lowerA, b); // calculate b_new
    }

    VectorXd
    calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd DInverse, VectorXd bNew, VectorXd x, double w) {
        return forward_subst(diagonalA + w*lowerA, (1-w)*diagonalA*x -  w* upperA*x) + bNew;
    }
};


VectorXd Gauss_Siedel_Iteration(const MatrixXd &A, const VectorXd &b, const VectorXd &x0, double tol)
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

//    std::cout << "Original:" <<std::endl << copiedA << std::endl;
//    std::cout << "Lower:" << std::endl << lower_A  << std::endl;
//    std::cout << "Upper" << std::endl << upper_A << std::endl;
//    std::cout << "Diagonal" << std::endl << diagonal_A << std::endl;
//    std::cout << lower_A + upper_A + diagonal_A - copiedA << std::endl;

    // Preparation for calculation
    MatrixXd D_inverse = diagonal_A.inverse();

    VectorXd b_new = forward_subst(diagonal_A + lower_A, b); // calculate b_new

    int ic = 0; // iteration count

    while (residual.norm() > stopIterResidual) {
        x = -1*forward_subst(diagonal_A + lower_A, upper_A*x) + b_new;
        residual = b - copiedA * x;
        ic += 1;
    }
    std::cout << "Iter count: " << ic << std::endl;

//   #TODO: Return tuple here with iter count
    return x;

}

VectorXd SOR_Iteration(const MatrixXd &A, const VectorXd &b, const VectorXd &x0, double tol, double w)
// w represents the omega or weight. Should be between 0 and 2
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

//    std::cout << "Original:" <<std::endl << copiedA << std::endl;
//    std::cout << "Lower:" << std::endl << lower_A  << std::endl;
//    std::cout << "Upper" << std::endl << upper_A << std::endl;
//    std::cout << "Diagonal" << std::endl << diagonal_A << std::endl;
//    std::cout << lower_A + upper_A + diagonal_A - copiedA << std::endl;

    // Preparation for calculation
    MatrixXd D_inverse = diagonal_A.inverse();

    VectorXd b_new = w*forward_subst(diagonal_A + w*lower_A, b); // calculate b_new

    int ic = 0; // iteration count

    while (residual.norm() > stopIterResidual) {
        x = forward_subst(diagonal_A + w*lower_A, (1-w)*diagonal_A*x -  w* upper_A*x) + b_new;
        residual = b - copiedA * x;
        ic += 1;
    }
    std::cout << "Iter count: " << ic << std::endl;

//   #TODO: Return tuple here with iter count
    return x;

}

#endif //CPPCODETEST_ITERATIVEMETHODS_HPP
