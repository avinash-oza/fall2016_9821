//
// Created by avi on 12/3/16.
//

#include "BinomialTrees.hpp"

#ifndef CPPCODETEST_TRINOMIALTREES_HPP
#define CPPCODETEST_TRINOMIALTREES_HPP


class EuropeanTrinomialTreePricer : public EuropeanBinomialTreePricer
{
public:
    using EuropeanBinomialTreePricer::EuropeanBinomialTreePricer;

    virtual TREE_RESULT calculateTree(long N)
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
                double riskNeutralDiscountedValue = calculateRiskNeutralDiscountedValue(i, j, deltaT, u, d, pu, pm, pd, optionPrices);
//                optionPrices[i] = std::max(riskNeutralDiscountedValue, K - S*std::pow(u, j - i)* std::pow(d, i)); // american
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

    virtual double calculateRiskNeutralDiscountedValue(int i, int j, double deltaT, double u ,double d, double pu, double pm, double pd, const VectorXd &optionPrices)
    {
        return std::exp(-r*deltaT)* (optionPrices(i) * pu + pm * optionPrices(i+1) + optionPrices(i + 2) * pd);
    }

};

class AmericanTrinomialTreePricer : public EuropeanTrinomialTreePricer
{
public:
    using EuropeanTrinomialTreePricer:: EuropeanTrinomialTreePricer;
    virtual double calculateRiskNeutralDiscountedValue(int i, int j, double deltaT, double u, double d, double pu, double pm, double pd, const VectorXd &optionPrices)
    {
        double europeanPrice = EuropeanTrinomialTreePricer::calculateRiskNeutralDiscountedValue(i, j, deltaT, u, d, pu, pm, pd, optionPrices);
        double intrinsicValue =  K - S*std::pow(u, j - i);

        return std::max(europeanPrice, intrinsicValue);
    }

};

#endif //CPPCODETEST_TRINOMIALTREES_HPP
