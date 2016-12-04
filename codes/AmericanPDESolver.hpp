//
// Created by avi on 11/28/16.
//

#ifndef CPPCODETEST_AMERICANPDESOLVER_HPP
#define CPPCODETEST_AMERICANPDESOLVER_HPP

#include "EuropeanPDESolver.hpp"
#include "uFunctions.hpp"

class AmericanPDESolver : public EuropeanPDESolver
{
public:
    /**
     *
     * @param gLeftFunc : Should be gAmericanLeft function unless another function is required
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
    AmericanPDESolver(uOptionFunction &gLeftFunc, uOptionFunction &gRightFunc, uOptionFunction &f, double t0,
                         double S0, double K, double T, double q, double r, double sigma, int M, double alphatemp) :
            EuropeanPDESolver(gLeftFunc, gRightFunc, f, t0, S0, K, T, q, r, sigma, M, alphatemp) {};

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
        return calculateVapprox(approximations) + (blackScholesPrice - europeanFDPrice);
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

/**
 * Given a mesh of approximations, this method returns a vector of optimal exercise Spot prices.
 * It first finds S_Nopt, S_Nopt + 1, and then computes the average
 * @param approximations : Approximations from an FD method
 * @return : Vector of SOptimal values foir early exercise
 */
    VectorXd findSEarlyExerciseSoptimal(MatrixXd &approximations)
    {
        // given the output of a method, this returns a 2 column vector.
        // The first column is U(N_opt), second column is U(N_opt+1)
        VectorXd earlyExerciseSpotValues(M + 1);
        double tol = std::pow(10, -14);

        // At the boundary/maturity (t=T), there is no difference between early exercise and holding
        earlyExerciseSpotValues(0) = K;

        // loop through each time vector and try to find the early exercise point
        for(int timeIndex = 1; timeIndex < M + 1; ++timeIndex)
        {
            VectorXd earlyExerciseValues = calculateEarlyExerciseVector(timeIndex);
            VectorXd uMeshValues = approximations.row(timeIndex);

            // loop through each x value in the time to see if it fills our criteria
            for(int xIndex = 0; xIndex < N - 1; ++xIndex)
            {
                // The logic does not make sense, but it seems to work
                if(uMeshValues(xIndex + 1) > earlyExerciseValues(xIndex))
                {
                    // this value is the early exercise value. Convert to spot price
                    double optimalValueLeft = K*std::exp(mesh.getX(xIndex));
                    double optimalValueRight = K*std::exp(mesh.getX(xIndex +1));
                    // find average
                    earlyExerciseSpotValues(timeIndex) = (optimalValueLeft + optimalValueRight)/ 2.0;
//                     to print the time in S,t and the early exercise value found
//                    std::cout << T - 2*mesh.getT(timeIndex)/(sigma*sigma) << "," << earlyExerciseSpotValues(timeIndex) << std::endl;
                    break;
                }
            }
        }
        return earlyExerciseSpotValues;
    }
};


#endif //CPPCODETEST_AMERICANPDESOLVER_HPP
