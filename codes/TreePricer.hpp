//
// Created by avi on 12/3/16.
//

#ifndef CPPCODETEST_TREEPRICER_HPP
#define CPPCODETEST_TREEPRICER_HPP

#include "BinomialTrees.hpp"
using namespace Eigen;

typedef std::tuple<double, double, double, double, double, Eigen::VectorXd> TREE_RESULT;
enum OPTION_TYPE{CALL, PUT};

class TreePricer

{

public:
    TreePricer(double S, double K, double T, double q, double r, double sigma, OPTION_TYPE optionType) : S(S), K(K), T(T), q(q), r(r),
                                                                                 sigma(sigma), optionType(optionType) {}


    MatrixXd calculateTreeforNInterations(long startN, long endN, double exactPrice)
    {
        long numberOfIterations = endN - startN;
        MatrixXd results = MatrixXd::Zero(numberOfIterations + 1, 3);

        for(long i = startN; i <= endN ;++i)
        {
            long matrixLocation = i - startN; // remap the values to start from the 0th location
            TREE_RESULT pricerResult = calculateTree(i);
            results(matrixLocation, 0) = i;
            results(matrixLocation, 1) = extractPrice(pricerResult);
            results(matrixLocation, 2) = calculateApproximationError(pricerResult, exactPrice);
            if(i % 20 == 0)
            {
                std::cout << "Priced " << i << std::endl;
            }
        }

        return results;
    }

    void calculateMinApproximationErrors(const MatrixXd &treeResults, double exactValue, long startN, long endN)
    {
        long numberOfIterations = endN - startN;
        std::cout << "MIN APPROX ERRORS:" << std::endl;
        for(long i = 1; i <= numberOfIterations - 1 ;++i)
        {
            double currentApproximationError = treeResults(i, 2);
            double nextApproximationError = treeResults(i + 1, 2);

            if(currentApproximationError < nextApproximationError)
            {
                std::cout << treeResults(i, 0)
                          << ","
                          << treeResults(i, 1)
                          << ","
                          << treeResults(i, 2)
                          << std::endl;
            }

        }
        std::cout << "END MIN APPROX ERRORS" << std::endl;
    }

    virtual TREE_RESULT calculateTree(long N)
    {
        // dummy function that is overridden
        return std::make_tuple(0.0, 0.0,0.0, 0.0, 0.0, VectorXd::Zero(10));
    }


    double calculateApproximationError(TREE_RESULT treeResult, double exactPrice)
    {
        double price = extractPrice(treeResult);
        return std::abs(price - exactPrice);
    }

    double extractPrice(TREE_RESULT binomialTreeResult)
    {
        return std::get<0>(binomialTreeResult);
    }

    double extractDelta(TREE_RESULT binomialTreeResult)
    {
        return std::get<1>(binomialTreeResult);
    }

    double extractGamma(TREE_RESULT binomialTreeResult)
    {
        return std::get<2>(binomialTreeResult);
    }

    double extractTheta(TREE_RESULT treeResult)
    {
        return std::get<3>(treeResult);
    }

    virtual long double
    calculateRiskNeutralDiscountedValue(double deltaT, double p, const VectorXd &optionPrices, int i, double u, int j) const {
        // calculates european value by default
        return exp(-r * deltaT) * (optionPrices[i] * p + optionPrices[i + 1] * (1 - p));
    }

    virtual long double calculateIntrinsicValue(int N, double u, int i, double d)
    {
        return 0.0; // for the european case there is no intrinsic value
    }

    virtual const double getFinalOptionPrice(double u, double d, int j, int i, double riskNeutralDiscountedValue) const {
        return riskNeutralDiscountedValue;
    }

    double calculatePayoff(double S) const
    {
        return optionType == CALL ? calculateCallPayoff(S) : calculatePutPayoff(S);
    }

    double varianceReductionPrice(TREE_RESULT americanTreeResult, TREE_RESULT  europeanTreeResult, BlackScholesOption &blackScholesOption)
    {
        return extractPrice(americanTreeResult) - extractPrice(europeanTreeResult) + blackScholesOption.price();
    }

    double varianceReductionDelta(TREE_RESULT americanTreeResult, TREE_RESULT  europeanTreeResult, BlackScholesOption &blackScholesOption)
    {
        return extractDelta(americanTreeResult) - extractDelta(europeanTreeResult) + blackScholesOption.delta();
    }

    double varianceReductionGamma(TREE_RESULT americanTreeResult, TREE_RESULT  europeanTreeResult, BlackScholesOption &blackScholesOption)
    {
        return extractGamma(americanTreeResult) - extractGamma(europeanTreeResult) + blackScholesOption.gamma();
    }

    double varianceReductionTheta(TREE_RESULT americanTreeResult, TREE_RESULT  europeanTreeResult, BlackScholesOption &blackScholesOption)
    {
        return extractTheta(americanTreeResult) - extractTheta(europeanTreeResult) + blackScholesOption.theta();
    }

    /**
     *
     * @param optimalN : optimal N found from somewhere else
     *@param x00: Initial guess left
     * @param x0 : initial guess right
     * @param consecutiveApproximation the tolerance between 2 approximations to consider it as converged
     * @return
     */
    double calculateImpliedVolatilityViaSecantMethod(long optimalN, double marketPrice, double x00, double x0,
                                                     double consecutiveApproximation)
    {
        double optimalNPrice = marketPrice;

        double xNew = x0;
        double xOld = x00;
        double xOldest = 0;
        int count = 1;

        double xOldVolatilityPrice = calculateTreeForNandSigma(optimalN, xOld );
        double xOldestVolatilityPrice = calculateTreeForNandSigma(optimalN, xOldest);
        std::cout << count << "\t" << xOld << std::endl;
        while(std::abs(xNew - xOld) > consecutiveApproximation)
        {
            ++count;
            xOldest = xOld;
            xOld = xNew;
            // update the prices for the next iteration
            xOldVolatilityPrice = calculateTreeForNandSigma(optimalN, xOld );
            xOldestVolatilityPrice = calculateTreeForNandSigma(optimalN, xOldest);
            // calculate next xNew
            xNew = xOld - (xOldVolatilityPrice - optimalNPrice) * (xOld -xOldest) / (xOldVolatilityPrice - xOldestVolatilityPrice);

            std::cout << count << "\t" << xOld << std::endl;
        }

        return xNew;
    }


protected:
    double S;
    double K;
    double T;
    double q;
    double r;
    double sigma;
    OPTION_TYPE optionType;

    // simple payoff functions
    double calculateCallPayoff(double S) const
    {
        return std::max(S - K, 0.0);
    }

    double calculatePutPayoff(double S) const
    {
        return std::max(K - S, 0.0);
    }

    /**
    * Helps to calculate option price for a sigma and N
    * @param N
    * @param sigma
    * @return double
    */
    double calculateTreeForNandSigma(long N, double _sigma)
    {
        double oldSigma = sigma; // store old sigma
        sigma = _sigma; // overwrite to ours
        TREE_RESULT treeResult = calculateTree(N);
        double price = extractPrice(treeResult);
        //restore back the old sigma
        sigma = oldSigma;
        return price;
    }

};


#endif //CPPCODETEST_TREEPRICER_HPP
