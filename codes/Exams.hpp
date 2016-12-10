//
// Created by Tushar Chawla on 12/9/16.
//

#ifndef CPPCODETEST_EXAMS_HPP
#define CPPCODETEST_EXAMS_HPP

#include <iostream>
#include "BinomialTrees.hpp"
#include "AmericanPDESolver.hpp"

namespace Exam2015

{
    namespace Question1 {


        double S = 47.0;
        double K = 50.0;
        double T = 7.0 / 12.0; // make sure to do double division
        double q = 0.01;
        double r = 0.03;
        double sigma = 0.25;
        OPTION_TYPE optionType = PUT;

        void Question1Part1() {


            AmericanBinomialTreePricer americanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            TREE_RESULT tree_result = americanBinomialTreePricer.calculateTree(9);
            std::cout << americanBinomialTreePricer.extractDelta(tree_result) << "\t";
            std::cout << americanBinomialTreePricer.extractGamma(tree_result) << "\t";
            std::cout << americanBinomialTreePricer.extractTheta(tree_result);

        }

        void Question1Part2() {


            AmericanBinomialTreePricer americanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            BlackScholesPutOption blackScholesPutOption(S, K, T, q, r, sigma);
            double previousValue = 0.0;

            int n = 0;

            for (int i = 1; i < 1000; ++i) {
                TREE_RESULT american_tree_result = americanBinomialTreePricer.calculateTree(i);
                TREE_RESULT european_tree_result = europeanBinomialTreePricer.calculateTree(i);

                double currentValue = americanBinomialTreePricer.varianceReductionPrice(american_tree_result,
                                                                                        european_tree_result,
                                                                                        blackScholesPutOption);
                if (std::abs(currentValue - previousValue) < 0.0001) break;

                previousValue = currentValue;
                std::cout << previousValue << std::endl;
                n = i;
            }
            std::cout << n << std::endl;

        }

        void Question1Part3() {
            EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            //TREE_RESULT european_tree_result = europeanBinomialTreePricer.calculateTree()
        }
    }

    namespace Question2
    {

        double S0 = 27.0;
        double K = 30.0;
        double T = 5.0 / 12.0; // make sure to do double division
        double q = 0.01;
        double r = 0.03;
        double sigma = 0.35;


        OPTION_TYPE optionType = PUT;
        long M = 4;
        double tol = std::pow(10, -6);
        double omega = 1.2;
        double alphatemp = 0.5;

            void Question2Part1()
            {
                fOption fOption(sigma, S0, q, K, T, r);
                gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
                gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
                gAmericanRight gRightOption(sigma, S0, q, K, T, r);
                gEuropeanRight gEuropeanRight(sigma, S0, q, K, T, r);


                for (int i = 1; i <= 4; ++i)
                {
                    M = pow(4,i);
                    AmericanPDESolver solverAmerican(gAmericanLeftFunc, gEuropeanRight, fOption, 0.0, S0, K, T, q, r,
                                                     sigma,
                                                     M, alphatemp);
                    solverAmerican.setUp();
                    std::cout << M << "\t" << alphatemp << "\t" << solverAmerican.getAlpha() << "\t "
                              << solverAmerican.getN() << "\t" <<
                              solverAmerican.get_xLeft() << "\t" << solverAmerican.get_xRight() << "\t"
                              << solverAmerican.getXCompute() <<  "\t" << solverAmerican.get_tFinal() <<
                              "\t" << solverAmerican.getDeltaTau() << "\t" << solverAmerican.getDeltaX() << std::endl;
                }
            }
        void Question2Part2()
        {

        }
        }
    }





#endif //CPPCODETEST_EXAMS_HPP
