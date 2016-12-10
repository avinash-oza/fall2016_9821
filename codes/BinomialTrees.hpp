//
// Created by avi on 12/3/16.
//

#ifndef CPPCODETEST_BINOMIALTREES_HPP
#define CPPCODETEST_BINOMIALTREES_HPP

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "BlackScholes.hpp"
#include "TreePricer.hpp"

enum BARRIER_TYPE {UP_AND_OUT, DOWN_AND_OUT};

using namespace Eigen;

class EuropeanBinomialTreePricer : public TreePricer
{
public:
    EuropeanBinomialTreePricer(double S, double K, double T, double q, double r, double sigma, OPTION_TYPE optionType) : TreePricer(S, K, T, q, r, sigma, optionType)  {};

    virtual TREE_RESULT calculateTree(long N)
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
//            optionPrices[i] = std::max(K - S* std::pow(u, N - i) * std::pow(d, i) , 0.0);
            optionPrices[i] = calculatePayoff(S* std::pow(u, N - i) * std::pow(d, i));
        }

        for (int j = N - 1; j >= 0; j--)
        {
            for(int i = 0; i < j + 1; i++)
            {
                double riskNeutralDiscountedValue = calculateRiskNeutralDiscountedValue(deltaT, p, optionPrices, i, u , j);
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



        //----- To print Tree

        /*std::cout << V22_P <<" " << V21_P << " " << V20_P << std::endl;
        std::cout << V11_P << " " << V10_P << " " << std::endl;
        std::cout << V00_P << std::endl;*/


        return std::make_tuple(optionPrices[0], Delta_P, Gamma_P, Theta_P, INT_MIN);
    }

    virtual  TREE_RESULT averageBinomialTree(long N)
        // N is the upper bound. Ex- input of 1280 will do 1280 and 1281
    {
        const TREE_RESULT binomialTreeNPlus1 = calculateTree(N + 1);
        const TREE_RESULT binomialTreeN = calculateTree(N);

        // take average of each result
        VectorXd result = VectorXd::Zero(4);

        result(0) = (std::get<0>(binomialTreeNPlus1) + std::get<0>(binomialTreeN)) / 2.0;
        result(1) = (std::get<1>(binomialTreeNPlus1) + std::get<1>(binomialTreeN)) / 2.0;
        result(2) = (std::get<2>(binomialTreeNPlus1) + std::get<2>(binomialTreeN)) / 2.0;
        result(3) = (std::get<3>(binomialTreeNPlus1) + std::get<3>(binomialTreeN)) / 2.0;

        return std::make_tuple(result(0), result(1), result(2), result(3), INT_MIN);
    }

    virtual TREE_RESULT calculateTreeBlackScholes(long N)
    {
        double deltaT = T/N;
        double u = std::exp(sigma* std::sqrt(deltaT));
        double d = 1.0/u;
        double p = (std::exp((r - q)*deltaT) - d)/(u - d);
        std::vector<double> optionPrices(N + 1);
        BlackScholesPutOption blackScholesOption(S, K, T, q, r, sigma);

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
            double intrinsicValue = calculateIntrinsicValue(N, u, i, d); // needed to calculate intrinsic value for american option
            double SPrice = S* std::pow(u, N - 1 - i) * std::pow(d, i);
            blackScholesOption.setS(SPrice);
            blackScholesOption.setT(deltaT);
            optionPrices[i] = std::max(intrinsicValue, blackScholesOption.price());

        }

        for (int j = N - 2; j >= 0; j--)
        {
            for(int i = 0; i < j + 1; i++)
            {
                double riskNeutralDiscountedValue = calculateRiskNeutralDiscountedValue(deltaT, p, optionPrices, i, u , j);
            optionPrices[i] = getFinalOptionPrice(u, d, j, i, riskNeutralDiscountedValue);

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
    }



    virtual TREE_RESULT BlackScholesWithRichardsonExtrapolation(long N)
    {
        auto bbs1 = calculateTreeBlackScholes(N);
        auto bbshalfN = calculateTreeBlackScholes(N / 2.0);

        // take average of each result
        VectorXd result = VectorXd::Zero(4);

        result(0) = 2.0*(std::get<0>(bbs1)) - std::get<0>(bbshalfN);
        result(1) = 2.0*(std::get<1>(bbs1)) - std::get<1>(bbshalfN);
        result(2) = 2.0*(std::get<2>(bbs1)) - std::get<2>(bbshalfN);
        result(3) = 2.0*(std::get<3>(bbs1)) - std::get<3>(bbshalfN);

        return std::make_tuple(result(0), result(1), result(2), result(3), INT_MIN);
    }


};

class AmericanBinomialTreePricer : public EuropeanBinomialTreePricer
{
public:
    AmericanBinomialTreePricer(double S, double K, double T, double q, double r, double sigma, OPTION_TYPE optionType)
            : EuropeanBinomialTreePricer(S, K, T, q, r, sigma, optionType) {}

    virtual long double
    calculateRiskNeutralDiscountedValue(double deltaT, double p, const std::vector<double> &optionPrices, int i, double u , int j) const {
        long double europeanPrice = TreePricer::calculateRiskNeutralDiscountedValue(deltaT, p, optionPrices, i, u , j);
//        long double intrinsicValue =  K - S*std::pow(u, j - i)* std::pow(1.0/u, i);
        long double intrinsicValue =  calculatePayoff(S*std::pow(u, j - i)* std::pow(1.0/u, i));

        return std::max(europeanPrice, intrinsicValue);
    }

    virtual long double calculateIntrinsicValue(int N, double u, int i, double d)
    {
        return calculatePayoff(S* std::pow(u, N - i - 1) * std::pow(d, i)); // for the european case there is no intrinsic value
    }

    virtual const double getFinalOptionPrice(double u, double d, int j, int i, double riskNeutralDiscountedValue) const {
        double payoff = calculatePayoff(S * pow(u, j - i) * pow(d, i));
        return std::max(riskNeutralDiscountedValue, payoff);
    }
};


class BarrierOptionBinomialTreePricer : public EuropeanBinomialTreePricer
{
public:
    BarrierOptionBinomialTreePricer(double S, double K, double T, double q, double r, double sigma, double B, BARRIER_TYPE barrierType, OPTION_TYPE optionType)
            : EuropeanBinomialTreePricer(S, K, T, q, r, sigma, optionType), B(B), barrierType(barrierType) {}

    virtual TREE_RESULT calculateTree(long N) {
        double deltaT = T / N;
        double u = std::exp(sigma * std::sqrt(deltaT));
        double d = 1 / u;
        double p = (std::exp((r - q) * deltaT) - d) / (u - d);
        std::vector<double> optionPrices(N + 1);

        double V22_P, V21_P, V20_P, V22_C, V21_C, V20_C, V11_P, V10_P, V11_C, V10_C, V00_P, V00_C;
        V22_P = V21_P = V20_P = V22_C = V21_C = V20_C = V11_P = V10_P = V11_C = V10_C = V00_P = V00_C = 0;
        double S11, S10, S22, S21, S20;
        double Delta_P, Delta_C, Gamma_P, Gamma_C, Theta_C, Theta_P;

        for (int i = 0; i < N; i++) {
            optionPrices[i] = 0.0;
        }

        double payoff;
        double currentSpot;
        for (int i = 0; i < N + 1; i++) {
            currentSpot = S * std::pow(u, N - i) * std::pow(d, i);
            double unHitPayoff = std::max(currentSpot - K, 0.0); //CHANGE FOR CALL OR PUT
            payoff = calculateFinalOptionValue(currentSpot, unHitPayoff); // checks to see if the barrier has been hit

            optionPrices[i] = payoff;
        }

        for (long j = N - 1; j >= 0; j--) {
            for (int i = 0; i < j + 1; i++) {
                currentSpot = S * std::pow(u, j - i) * std::pow(d, i);

                double riskNeutralDiscountedValue = calculateRiskNeutralDiscountedValue(deltaT, p, optionPrices, i, u, j);
                optionPrices[i] = calculateFinalOptionValue(currentSpot, riskNeutralDiscountedValue); // european

                if (j == 2) {
                    V22_P = optionPrices[0];
                    V21_P = optionPrices[1];
                    V20_P = optionPrices[2];

//                V22_C = v_call(0);
//                V21_C = v_call(1);
//                V20_C = v_call(2);
                } else if (j == 1) {
                    V11_P = optionPrices[0];
                    V10_P = optionPrices[1];

//                V11_C = v_call(0);
//                V10_C = v_call(1);
                } else if (j == 0) {
                    V00_P = optionPrices[0];

//                V00_C = v_call(1);
                }
            }

        }
        S11 = u * S;
        S10 = d * S;
        S22 = u * u * S;
        S21 = u * d * S;
        S20 = d * d * S;

        Delta_P = (V10_P - V11_P) / (S10 - S11);
        //        Delta_C = (V10_C - V11_C) / (S10 - S11);
        Gamma_P = ((V20_P - V21_P) / (S20 - S21) - (V21_P - V22_P) / (S21 - S22)) / ((S20 - S22) / 2.0);
        //        Gamma_C = ((V20_C - V21_C) / (S20 - S21) - (V21_C - V22_C) / (S21 - S22)) / ((S20 - S22) / 2);
        Theta_P = (V21_P - V00_P) / (2.0 * deltaT);
        //        Theta_C = (V21_C - V00_C) / (2 * deltaT);


        return std::make_tuple(optionPrices[0], Delta_P, Gamma_P, Theta_P, INT_MIN);
    }

    double calculateFinalOptionValue(double currentSpot, double unHitPayoff) const {
        // unHitPayoff - If the barrier was not hit, return this value
        double payoff;
        if(barrierType == DOWN_AND_OUT && (currentSpot < B))
            {
                //option knocked out
                payoff = 0.0;

            }
            else if(barrierType == UP_AND_OUT && (currentSpot > B))
            {
                //option knocked out
                payoff = 0.0;

            }
            else{
                payoff = unHitPayoff;
            }
        return payoff;
    }

protected:
    double B;
    BARRIER_TYPE barrierType;
};

#endif //CPPCODETEST_BINOMIALTREES_HPP
