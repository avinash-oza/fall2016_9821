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
        virtual VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd Dinverse, VectorXd b, VectorXd x0, double w)=0;
        virtual VectorXd calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd Dinverse, VectorXd bNew, VectorXd x, double w)=0;
};

class JacobiIteration: public IterationMethod
{
public:
    VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd D_inverse, VectorXd b, VectorXd x0, double w) {
        return D_inverse*b;
    }

    VectorXd calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd DInverse, VectorXd bNew, VectorXd x, double w){
        return -1* DInverse * (lowerA * x + upperA * x) + bNew;
    }
};

class GaussSiedelIteration: public IterationMethod
{
public:
    VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd Dinverse, VectorXd b, VectorXd x0, double w) {
        return forward_subst(diagonalA + lowerA, b); // calculate b_new;
    }

    VectorXd calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd DInverse, VectorXd bNew, VectorXd x, double w) {
        return -1*forward_subst(diagonalA + lowerA, upperA*x) + bNew;
    }
};


class SORIteration: public IterationMethod
{
public:
    VectorXd calculateBnew(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd D_inverse, VectorXd b, VectorXd x0, double w) {
//        return w*forward_subst(diagonalA + w*lowerA, b); // calculate b_new
        return b;
    }

    virtual double calculateSORValue(double currentValue, int currentIndex)
    {
        // for European options we don't need to edit the value
        return currentValue;
    }

    virtual VectorXd calculateIteration(MatrixXd lowerA, MatrixXd upperA, MatrixXd diagonalA, MatrixXd DInverse, VectorXd bNew, VectorXd x, double w) {
        MatrixXd A = lowerA + upperA + diagonalA;
        long size = x.rows();
        VectorXd xSOR = VectorXd::Zero(size);
        double currentValue;

        for( int j = 0; j < size; ++j)
        {
            if (j == 0)
            {
                currentValue = (1 - w) * x(j) + w / A(j,j) * (bNew(j) - A(j,j+1) * x(j+1));
            }
            else
            {
                if (j == (size - 1))
                {
                    currentValue = (1 - w) * x(j) + w / A(j,j) * (bNew(j) - A(j,j-1) * xSOR(j-1));
                }
                else
                {
                    currentValue = (1 - w) * x(j) + w / A(j,j) * (bNew(j) - A(j,j-1) * xSOR(j-1) - A(j,j+1) * x(j+1));
                }
            }
            xSOR(j) = calculateSORValue(currentValue, j);

        }
        return xSOR;
    }

};

class ProjectedSORIteration : public SORIteration
{
public:
    ProjectedSORIteration(VectorXd & _earlyExerciseValues) : earlyExerciseValues(_earlyExerciseValues) {};

    virtual double calculateSORValue(double currentValue, int currentIndex) override {
        double &earlyExerciseValue = earlyExerciseValues(currentIndex);
        return std::max(earlyExerciseValue, currentValue);;
    }

private:
    VectorXd earlyExerciseValues;
};

#endif //CPPCODETEST_ITERATIVEMETHODS_HPP
