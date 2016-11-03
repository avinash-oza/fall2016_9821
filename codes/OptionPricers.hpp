//
// Created by avi on 9/25/16.
//

#ifndef CPPCODETEST_OPTIONPRICERS_HPP
#define CPPCODETEST_OPTIONPRICERS_HPP
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "BlackScholes.hpp"


std::tuple<double, double, double, double, double> binomialTree(double S, double K, double T, double q, double r, int N,
                                                        double sigma)
//S is spot price, K is strike price, T is time to maturity in years, q is the continuous dividend, r is the risk free rate
// N is the number of time steps
{
    double deltaT = T/N;
    double u = std::exp(sigma* std::sqrt(deltaT));
    double d = 1/u;
    double p = (std::exp((r - q)*deltaT) - d)/(u - d);
    std::vector<double> optionPrices(N + 1);

    double V22_P, V21_P, V20_P, V22_C, V21_C, V20_C, V11_P, V10_P, V11_C, V10_C, V00_P, V00_C;
    V22_P = V21_P = V20_P = V22_C = V21_C = V20_C = V11_P = V10_P = V11_C = V10_C = V00_P = V00_C = 0;
    double S11, S10, S22, S21, S20;
    double Delta_P, Delta_C, Gamma_P, Gamma_C, Theta_C, Theta_P;

    for (int i = 0; i< N; i++)
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
//            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i)); // american
            optionPrices[i] = riskNeutralDiscountedValue; // european

            if (j == 2)
            {
                V22_P = optionPrices[0];
                V21_P = optionPrices[1];
                V20_P = optionPrices[2];

//                V22_C = v_call(0);
//                V21_C = v_call(1);
//                V20_C = v_call(2);
            }
            else if (j == 1)
            {
                V11_P = optionPrices[0];
                V10_P = optionPrices[1];

//                V11_C = v_call(0);
//                V10_C = v_call(1);
            }
            else if (j == 0)
            {
                V00_P = optionPrices[0];

//                V00_C = v_call(1);
            }
        }

        }
    S11 = u*S;
    S10 = d*S;
    S22 = u*u*S;
    S21 = u*d*S;
    S20 = d*d*S;

    Delta_P = (V10_P - V11_P) / (S10 - S11);
    //        Delta_C = (V10_C - V11_C) / (S10 - S11);
    Gamma_P = ((V20_P - V21_P) / (S20 - S21) - (V21_P - V22_P) / (S21 - S22)) / ((S20 - S22) / 2.0);
    //        Gamma_C = ((V20_C - V21_C) / (S20 - S21) - (V21_C - V22_C) / (S21 - S22)) / ((S20 - S22) / 2);
    Theta_P = (V21_P - V00_P) / (2.0 * deltaT);
    //        Theta_C = (V21_C - V00_C) / (2 * deltaT);

    return std::make_tuple(optionPrices[0], Delta_P, Gamma_P, Theta_P, INT_MIN);
    }

std::tuple<double, double, double, double, double> binomialBlackScholes(double S, double K, double T, double q, double r, int N,
                                                                double sigma)
//S is spot price, K is strike price, T is time to maturity in years, q is the continuous dividend, r is the risk free rate
// N is the number of time steps
{
    double deltaT = T/N;
    double u = std::exp(sigma* std::sqrt(deltaT));
    double d = 1.0/u;
    double p = (std::exp((r - q)*deltaT) - d)/(u - d);
    std::vector<double> optionPrices(N + 1);
    BlackScholesOption blackScholesOption;

    double V22_P, V21_P, V20_P, V22_C, V21_C, V20_C, V11_P, V10_P, V11_C, V10_C, V00_P, V00_C;
    V22_P = V21_P = V20_P = V22_C = V21_C = V20_C = V11_P = V10_P = V11_C = V10_C = V00_P = V00_C = 0;
    double S11, S10, S22, S21, S20;
    double Delta_P, Delta_C, Gamma_P, Gamma_C, Theta_C, Theta_P;

    for (int i = 0; i< N; i++)
    {
        optionPrices[i] = 0.0;
    }

    for (int i = 0 ; i < N; i++)
    {
//        optionPrices[i] = std::max(K - S* std::pow(u, N - i) * std::pow(d, i) , 0.0)
        double SPrice = S* std::pow(u, N - 1 - i) * std::pow(d, i);
        optionPrices[i] = blackScholesOption.putPrice(SPrice, K, deltaT, q, r, N, sigma);
    }

    for (int j = N - 2; j >= 0; j--)
    {
        for(int i = 0; i < j + 1; i++)
        {
            double riskNeutralDiscountedValue = std::exp(-r*deltaT)* (optionPrices[i] * p + optionPrices[i + 1] * (1 - p));
//            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i)); // american
            optionPrices[i] = riskNeutralDiscountedValue; // european

            if (j == 2)
            {
                V22_P = optionPrices[0];
                V21_P = optionPrices[1];
                V20_P = optionPrices[2];

//                V22_C = v_call(0);
//                V21_C = v_call(1);
//                V20_C = v_call(2);
            }
            else if (j == 1)
            {
                V11_P = optionPrices[0];
                V10_P = optionPrices[1];

//                V11_C = v_call(0);
//                V10_C = v_call(1);
            }
            else if (j == 0)
            {
                V00_P = optionPrices[0];

//                V00_C = v_call(1);
            }

        }
    }

    S11 = u*S;
    S10 = d*S;
    S22 = u*u*S;
    S21 = u*d*S;
    S20 = d*d*S;

    Delta_P = (V10_P - V11_P) / (S10 - S11);
    //        Delta_C = (V10_C - V11_C) / (S10 - S11);
    Gamma_P = ((V20_P - V21_P) / (S20 - S21) - (V21_P - V22_P) / (S21 - S22)) / ((S20 - S22) / 2);
    //        Gamma_C = ((V20_C - V21_C) / (S20 - S21) - (V21_C - V22_C) / (S21 - S22)) / ((S20 - S22) / 2);
    Theta_P = (V21_P - V00_P) / (2 * deltaT);
    //        Theta_C = (V21_C - V00_C) / (2 * deltaT);
    return std::make_tuple(optionPrices[0], Delta_P, Gamma_P, Theta_P, INT_MIN);
//    return optionPrices[0];
}

double averageBinomialTree(double S, double K, double T, double q, double r, int N, double sigma)
// N is the upper bound. Ex- input of 1280 will do 1280 and 1281
{
    const std::tuple<double, double, double, double, double> binomialTreeNPlus1 = binomialTree(S, K, T, q, r, N + 1, sigma);
    const std::tuple<double, double, double, double, double> binomialTreeN = binomialTree(S, K, T, q, r, N , sigma);


    double v1 = std::get<3>(binomialTreeNPlus1);
    double v2 = std::get<3>(binomialTreeN);
    return (v1 + v2) / 2.0;
}

double binomialBlackScholeswithRichardsonExtrapolation(double S, double K, double T, double q, double r, int N, double sigma)
{
    auto bbs1 = binomialBlackScholes(S,K, T, q, r, N, sigma);
    auto bbshalfN = binomialBlackScholes(S,K, T, q, r, N/2.0, sigma);

    return 2.0*(std::get<3>(bbs1)) - std::get<3>(bbshalfN);
}

void calculateTreesForN(int N)
{
    double S = 41.0;
    double K = 40.0;
    double T = 1.0;
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;
    BlackScholesOption blackScholesOption;

    std::tuple<double, double, double, double, double> binomialTreePrice = binomialTree(S, K, T, q, r, N, sigma);
    double averageBinomialTreePrice = averageBinomialTree(S, K, T, q, r, N, sigma);
    double blackScholesValue = blackScholesOption.putPrice(S, K, T, q, r, N, sigma);
    double blackScholeswithRichardsonExtrapolation = binomialBlackScholeswithRichardsonExtrapolation(S, K, T, q, r, N, sigma);
    double binomialBlackScholesPrice = std::get<3>(binomialBlackScholes(S, K, T, q, r, N, sigma));
    double richardsonPrice = binomialBlackScholeswithRichardsonExtrapolation(S, K, T, q, r, N, sigma);

//    double optionPrice = std::get<0>(binomialTreePrice);
    double optionPrice = richardsonPrice;
    double absDiff = std::abs(richardsonPrice - blackScholesValue);
    double deltaPut = blackScholesOption.putDelta(S, K, T, q, r, N, sigma);
//    std::cout << blackScholesValue << std::endl;
    std::cout << std::setprecision(12) << deltaPut << ","
//              << absDiff << ","
//              << N*absDiff << ","
//              << N*N*absDiff
              << std::endl;

}

std::tuple<double, double, double, double, double> trinomialTree(double S, double K, double T, double q, double r, int N,
                                                        double sigma)
//S is spot price, K is strike price, T is time to maturity in years, q is the continuous dividend, r is the risk free rate
// N is the number of time steps
{
    double deltaT = (T*1.0)/N;
    double u = std::exp(sigma*std::sqrt(3.0*deltaT));
    double d = 1.0/u;
    double pm = 2.0/3.0;
    double pu = 1.0/6.0 + (r - q - sigma*sigma/2)*std::sqrt(deltaT/(12.0*sigma*sigma));
    double pd = 1.0/6.0 - (r - q - sigma*sigma/2)*std::sqrt(deltaT/(12.0*sigma*sigma));

    double S_10, S_11, S_12, S_20, S_22, S_24;
    S_10 = S_11 = S_12 = S_20 = S_22 = S_24 = 0;


    MatrixXd importantValues(3, 5);
    importantValues.setZero();
    double Delta_P, Delta_C, Gamma_P, Gamma_C, Theta_C, Theta_P;

    Eigen::VectorXd optionPrices(2*N + 1);
    optionPrices.setZero();

    for (int i = 0 ; i < 2*N + 1; i++)
    {
        optionPrices(i) = std::max(K - S* std::pow(u, N - i) , 0.0);
    }

    for (int j = N - 1; j >= 0; j--)
    {
        for(int i = 0; i < 2*j + 1; i++)
        {
            double riskNeutralDiscountedValue = std::exp(-r*deltaT)* (optionPrices(i) * pu + pm * optionPrices(i+1) + optionPrices(i + 2) * pd);
//            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i)); // american
            optionPrices(i) = riskNeutralDiscountedValue; // european

            if (j == 2)
            {
                importantValues(2, 4) = optionPrices(0);
                importantValues(2, 2) = optionPrices(2);
                importantValues(2, 0) = optionPrices(4);

            }
            else if (j == 1)
            {
                importantValues(1, 2) = optionPrices(0);
                importantValues(1, 1) = optionPrices(1);
                importantValues(1, 0) = optionPrices(2);
            }
            else if (j == 0)
            {
                importantValues(0, 0) = optionPrices(0);

//                V00_C = v_call(1);
            }
        }

        }

    S_10 = d*S;
    S_12 = u*S;
    S_24 = u*u*S;
    S_22 = u*d*S;
    S_20 = d*d*S;

    Delta_P = (importantValues(1, 0) - importantValues(1, 2)) / (S_10 - S_12);
    //        Delta_C = (V10_C - V11_C) / (S10 - S11);
    Gamma_P = ((importantValues(2, 0) - importantValues(2, 2))/ (S_20 - S_22) - (importantValues(2, 2) - importantValues(2,4))/(S_22 - S_24))/ (S_10 - S_12);
    //        Gamma_C = ((V20_C - V21_C) / (S20 - S21) - (V21_C - V22_C) / (S21 - S22)) / ((S20 - S22) / 2);
    Theta_P = (importantValues(1, 1) - importantValues(0, 0)) / deltaT;
    //        Theta_C = (V21_C - V00_C) / (2 * deltaT);

    return std::make_tuple(optionPrices(0), Delta_P, Gamma_P, Theta_P, INT_MIN);
    }

std::tuple<double, double, double, double, double> trinomialBlackScholes(double S, double K, double T, double q, double r, int N,
                                                                 double sigma)
//S is spot price, K is strike price, T is time to maturity in years, q is the continuous dividend, r is the risk free rate
// N is the number of time steps
{
    double deltaT = (T*1.0)/N;
    double u = std::exp(sigma*std::sqrt(3.0*deltaT));
    double d = 1.0/u;
    double pm = 2.0/3.0;
    double pu = 1.0/6.0 + (r - q - sigma*sigma/2)*std::sqrt(deltaT/(12.0*sigma*sigma));
    double pd = 1.0/6.0 - (r - q - sigma*sigma/2)*std::sqrt(deltaT/(12.0*sigma*sigma));
    BlackScholesOption blackScholesOption;

    double S_10, S_11, S_12, S_20, S_22, S_24;
    S_10 = S_11 = S_12 = S_20 = S_22 = S_24 = 0;


    MatrixXd importantValues(3, 5);
    importantValues.setZero();
    double Delta_P, Delta_C, Gamma_P, Gamma_C, Theta_C, Theta_P;

    Eigen::VectorXd optionPrices(2*N + 1);
    optionPrices.setZero();

    for (int i = 0 ; i < 2*N - 1; i++)
    {
        double SPrice = S* std::pow(u, N - 1 - i) ;
        optionPrices(i) = blackScholesOption.putPrice(SPrice, K, deltaT, q, r, N, sigma);
    }

    for (int j = N - 2; j >= 0; j--)
    {
        for(int i = 0; i < 2*j + 1; i++)
        {
            double riskNeutralDiscountedValue = std::exp(-r*deltaT)* (optionPrices(i) * pu + pm * optionPrices(i+1) + optionPrices(i + 2) * pd);
//            optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i)); // american
            optionPrices(i) = riskNeutralDiscountedValue; // european

            if (j == 2)
            {
                importantValues(2, 4) = optionPrices(0);
                importantValues(2, 2) = optionPrices(2);
                importantValues(2, 0) = optionPrices(4);

            }
            else if (j == 1)
            {
                importantValues(1, 2) = optionPrices(0);
                importantValues(1, 1) = optionPrices(1);
                importantValues(1, 0) = optionPrices(2);
            }
            else if (j == 0)
            {
                importantValues(0, 0) = optionPrices(0);

//                V00_C = v_call(1);
            }
        }

    }

    S_10 = d*S;
    S_12 = u*S;
    S_24 = u*u*S;
    S_22 = u*d*S;
    S_20 = d*d*S;

    Delta_P = (importantValues(1, 0) - importantValues(1, 2)) / (S_10 - S_12);
    //        Delta_C = (V10_C - V11_C) / (S10 - S11);
    Gamma_P = ((importantValues(2, 0) - importantValues(2, 2))/ (S_20 - S_22) - (importantValues(2, 2) - importantValues(2,4))/(S_22 - S_24))/ (S_10 - S_12);
    //        Gamma_C = ((V20_C - V21_C) / (S20 - S21) - (V21_C - V22_C) / (S21 - S22)) / ((S20 - S22) / 2);
    Theta_P = (importantValues(1, 1) - importantValues(0, 0)) / deltaT;
    //        Theta_C = (V21_C - V00_C) / (2 * deltaT);

    return std::make_tuple(optionPrices(0), Delta_P, Gamma_P, Theta_P, INT_MIN);
}

std::tuple<double, double, double, double, double> trinomialBlackScholeswithRichardsonExtrapolation(double S, double K, double T, double q, double r, int N, double sigma)
{
    auto bbs1 = trinomialBlackScholes(S,K, T, q, r, N, sigma);
    auto bbshalfN = trinomialBlackScholes(S,K, T, q, r, N/2.0, sigma);
    double price = 2.0*(std::get<0>(bbs1)) - std::get<0>(bbshalfN);
    double delta = 2.0*(std::get<1>(bbs1)) - std::get<1>(bbshalfN);
    double gamma = 2.0*(std::get<2>(bbs1)) - std::get<2>(bbshalfN);
    double theta = 2.0*(std::get<3>(bbs1)) - std::get<3>(bbshalfN);

    return std::make_tuple(price, delta, gamma, theta, INT_MIN);
}


void calculateTrinomialTreesForN(int N)
{
    double S = 41.0;
    double K = 40.0;
    double T = 1.0;
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;
    BlackScholesOption blackScholesOption;

    std::tuple<double, double, double, double, double> binomialTreePrice = trinomialTree(S, K, T, q, r, N, sigma);
    double trinomialTreePrice = std::get<0>(binomialTreePrice);
    double deltaTrinomial = std::get<1>(binomialTreePrice);
    double gammaTrinomial = std::get<2>(binomialTreePrice);
    double thetaTrinomial = std::get<3>(binomialTreePrice);
    double blackScholesValue = blackScholesOption.putPrice(S, K, T, q, r, N, sigma);

//    double optionPrice = std::get<0>(binomialTreePrice);
//    double optionPrice = richardsonPrice;
//    double absDiff = std::abs(richardsonPrice - blackScholesValue);
//    double deltaPut = deltaBSPut(S, K, T, q, r, N, sigma);
//    std::cout << blackScholesValue << std::endl;
    std::cout << std::setprecision(12) << trinomialTreePrice << ","
              << deltaTrinomial << ","
            << gammaTrinomial << ","
              << thetaTrinomial << ","
//              << absDiff << ","
//              << N*absDiff << ","
//              << N*N*absDiff
              << std::endl;



}



#endif //CPPCODETEST_OPTIONPRICERS_HPP
