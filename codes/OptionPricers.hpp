//
// Created by avi on 9/25/16.
//

#ifndef CPPCODETEST_OPTIONPRICERS_HPP
#define CPPCODETEST_OPTIONPRICERS_HPP
#include <vector>
#include <cmath>


double binomialTree(double S, double K, double T, double q, double r, int N, double sigma)
//S is spot price, K is strike price, T is time to maturity in years, q is the continuous dividend, r is the risk free rate
// N is the number of time steps
{
    double deltaT = T/N;
    double u = std::exp(sigma* std::sqrt(deltaT));
    double d = 1/u;
    double p = (std::exp((r - q)*deltaT) - d)/(u - d);
    std::vector<double> optionPrices(N);
    for (int i = 1; i<N; i++)
    {
        optionPrices[i] = 0.0;
    }

    for (int i = 0 ; i < N; i++)
    {
        optionPrices[i] = std::max(K - S* std::pow(u, N - i) * std::pow(d, i) , 0.0);
    }
    for (int j = N ; j > 0; j--)
    {
        for(volatile int i = 0; i < j; i++)
        {
            double riskNeutralDiscountedValue = std::exp(-r*deltaT)* (optionPrices[i] * p + optionPrices[i + 1] * (1 - p));
            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i));
        }
    }
    return optionPrices[0];
}

#endif //CPPCODETEST_OPTIONPRICERS_HPP
