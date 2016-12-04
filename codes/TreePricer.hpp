//
// Created by avi on 12/3/16.
//

#ifndef CPPCODETEST_TREEPRICER_HPP
#define CPPCODETEST_TREEPRICER_HPP
typedef std::tuple<double, double, double, double, double> TREE_RESULT;

#include "BinomialTrees.hpp"


class TreePricer

{

public:
    TreePricer(double S, double K, double T, double q, double r, double sigma) : S(S), K(K), T(T), q(q), r(r),
                                                                                 sigma(sigma) {}

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
    calculateRiskNeutralDiscountedValue(double deltaT, double p, const std::vector<double> &optionPrices, int i, double u, int j) const {
        // calculates european value by default
        return exp(-r * deltaT) * (optionPrices[i] * p + optionPrices[i + 1] * (1 - p));
    }
protected:
    double S;
    double K;
    double T;
    double q;
    double r;
    double sigma;
};


#endif //CPPCODETEST_TREEPRICER_HPP
