//
// Created by avi on 9/25/16.
//

#ifndef CPPCODETEST_OPTIONPRICERS_HPP
#define CPPCODETEST_OPTIONPRICERS_HPP
#include <vector>
#include <cmath>
#include "BlackScholes.hpp"


double binomialTree(double S, double K, double T, double q, double r, int N, double sigma)
//S is spot price, K is strike price, T is time to maturity in years, q is the continuous dividend, r is the risk free rate
// N is the number of time steps
{
    double deltaT = T/N;
    double u = std::exp(sigma* std::sqrt(deltaT));
    double d = 1/u;
    double p = (std::exp((r - q)*deltaT) - d)/(u - d);
    std::vector<double> optionPrices(N + 1);
    for (int i = 1; i< N + 1; i++)
    {
        optionPrices[i] = 0.0;
    }

    for (int i = 0 ; i < N + 1; i++)
    {
        optionPrices[i] = std::max(K - S* std::pow(u, N - i) * std::pow(d, i) , 0.0);
    }

    for (int j = N - 1; j >= 0; j--)
    {
        for(int i = 0; i < j + 1; i++)
        {
            double riskNeutralDiscountedValue = std::exp(-r*deltaT)* (optionPrices[i] * p + optionPrices[i + 1] * (1 - p));
            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i)); // american
//            optionPrices[i] = riskNeutralDiscountedValue; // european
        }
    }
    return optionPrices[0];
}

double binomialBlackScholes(double S, double K, double T, double q, double r, int N, double sigma)
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
        double spotPrice = S*std::pow(u, N - 1 - i)* std::pow(d, i);
        optionPrices[i] = blackScholesPut(spotPrice, K, T - i*deltaT, q, r, N, sigma);
        // The put price
    }

    for (int j = N; j > 0; j--)
    {
        for(int i = 0; i < j; i++)
        {
            double riskNeutralDiscountedValue = std::exp(-r*deltaT)* (optionPrices[i] * p + optionPrices[i + 1] * (1 - p));
//            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i));
            optionPrices[i] = riskNeutralDiscountedValue;
        }
    }
    return optionPrices[0];
}

double averageBinomialTree(double S, double K, double T, double q, double r, int N, double sigma)
// N is the upper bound. Ex- input of 1280 will do 1280 and 1281
{
    return (binomialTree(S, K, T, q, r, N + 1, sigma) + binomialTree(S, K, T, q, r, N , sigma))/2.0;
}

double binomialBlackScholeswithRichardsonExtrapolation(double S, double K, double T, double q, double r, int N, double sigma)
{
    double bbs1 = blackScholesPut(S,K, T, q, r, N, sigma);
    double bbshalfN = blackScholesPut(S,K, T, q, r, N/2.0, sigma);

    return 2.0*bbs1 - bbshalfN;
}

void calculateTreesForN(int N)
{
    double S = 41.0;
    double K = 40.0;
    double T = 1;
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;

    double binomialTreePrice = binomialTree(S, K, T, q, r, N, sigma);
    double averageBinomialTreePrice = averageBinomialTree(S, K, T, q, r, N, sigma);
    double blackScholesValue = blackScholesPut(S, K, T, q, r, N, sigma);
    double blackScholeswithRichardsonExtrapolation = binomialBlackScholeswithRichardsonExtrapolation(S, K, T, q, r, N, sigma);
    double binomialBlackScholesPrice = binomialBlackScholes(S, K, T, q, r, N, sigma);

    double absDiff = std::abs(binomialTreePrice - blackScholesValue);
    std::cout << binomialTreePrice << "," << absDiff << "," << N*absDiff << "," << N*N*absDiff << std::endl;

}

#endif //CPPCODETEST_OPTIONPRICERS_HPP
