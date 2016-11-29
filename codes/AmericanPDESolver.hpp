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
};


#endif //CPPCODETEST_AMERICANPDESOLVER_HPP