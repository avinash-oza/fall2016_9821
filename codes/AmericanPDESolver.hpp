//
// Created by avi on 11/28/16.
//

#ifndef CPPCODETEST_AMERICANPDESOLVER_HPP
#define CPPCODETEST_AMERICANPDESOLVER_HPP

#include "EuropeanPDESolver.hpp"
#include "uFunctions.hpp"

class AmericanPutPDESolver : public EuropeanPutPDESolver
{
public:
    /**
     *
     * @param gLeftFunc : gLeft is hardcoded to the gAmericanLeft function
     * @param gRightFunc : right boundary function
     * @param f : boundary at the bottom of the rectangle
     * @param t0 : typically 0
     * @param S0 : inital spot price
     * @param K  : strike
     * @param T : maturity
     * @param q : dividend
     * @param r : risk free rate
     * @param sigma : volatility
     * @param M : number of intervals in x
     * @param alphatemp 
     */
    AmericanPutPDESolver(uOptionFunction &gLeftFunc, uOptionFunction &gRightFunc, uOptionFunction &f, double t0,
                         double S0, double K, double T, double q, double r, double sigma, int M, double alphatemp) :
            EuropeanPutPDESolver(gLeftFunc, gRightFunc, f, t0, S0, K, T, q, r, sigma, M, alphatemp)
    {

        // for the american option we have a constant left boundary:
        ///////THIS NEEDS TO BE CHECKED!!!!!!!!!!!!!!
        _gLeftFunc = gAmericanLeftFunc(sigma, S0, q, K, T, r);

    };

    double earlyExercisePremium(double x, double t)
    {
        return K*std::exp(a*x + b*t) * std::max(1.0 - std::exp(x), 0.0);
    }

    virtual double calculateUOnMesh(long timeIndex, long currentIndex,  double europeanUValue)
    {
        // here we need to take the max between the european and american value
        double x = mesh.getX(currentIndex);
        double t = mesh.getT(timeIndex);
        double earlyExerciseValue = earlyExercisePremium(x, t);

        return std::max(europeanUValue, earlyExerciseValue);
    }

    double priceVarianceReduction(MatrixXd &approximations, double europeanFDPrice, double blackScholesPrice)
    {
        return calculateVapprox(S0, approximations) + (blackScholesPrice - europeanFDPrice);
    }

    double calculateErrorPointWiseVarianceReduction(MatrixXd &approximations, double europeanFDPrice, double blackScholesPrice, double P_american_bin)
    {
        return std::abs(priceVarianceReduction(approximations, europeanFDPrice, blackScholesPrice) - P_american_bin);
    }

    VectorXd calculateEarlyExerciseVector(long timeIndex) {
        VectorXd earlyExerciseValues = VectorXd::Zero(N - 1);
        double x;
        double t = mesh.getT(timeIndex);

        for (int i = 1; i < N - 1; ++i) {
            x = mesh.getX(i);
            earlyExerciseValues(i - 1) = earlyExercisePremium(x, t);
        }
        return earlyExerciseValues;
    }

    virtual std::tuple<VectorXd, int>
    getSORIterationValues(long timeIndex, double tol, double omega, const VectorXd &U, const MatrixXd &A, const MatrixXd &B,
                          const MatrixXd &b) {
        // we need to use the projected SOR method in the American put calculation
        VectorXd earlyExerciseValues = calculateEarlyExerciseVector(timeIndex);
        ProjectedSORIteration iterationMethod(earlyExerciseValues);
        return consecutive_approximation_solver(A, B * U + b, earlyExerciseValues, tol, iterationMethod, omega);;
    }
};


#endif //CPPCODETEST_AMERICANPDESOLVER_HPP
