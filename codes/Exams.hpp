//
// Created by Tushar Chawla on 12/9/16.
// for

#ifndef CPPCODETEST_EXAMS_HPP
#define CPPCODETEST_EXAMS_HPP

#include <iostream>
#include "BinomialTrees.hpp"
#include "AmericanPDESolver.hpp"
#include "BarrierOption.hpp"
#include "BarrierPDE.hpp"

namespace Exam2015 {
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

//                double currentValue = americanBinomialTreePricer.varianceReductionPrice(american_tree_result,
//                                                                                        european_tree_result,
//                                                                                        blackScholesPutOption);
                double currentValue = europeanBinomialTreePricer.extractPrice(american_tree_result);
                if (std::abs(currentValue - previousValue) < 0.0001)
                {
                    n += 1; // increment by 1 as the optimial iteration is the one greater than this
                    break;
                }

                previousValue = currentValue;
//                std::cout << previousValue << std::endl;
                n = i;
            }
            std::cout << n << std::endl;

        }

        void Question1Part3() {
            EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            //TREE_RESULT european_tree_result = europeanBinomialTreePricer.calculateTree()
        }
    }

    namespace Question2 {

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

        void Question2Part1() {
            fOption fOption(sigma, S0, q, K, T, r);
            gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
            gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
            gAmericanRight gRightOption(sigma, S0, q, K, T, r);
            gEuropeanRight gEuropeanRight(sigma, S0, q, K, T, r);


            for (int i = 1; i <= 4; ++i) {

                AmericanPDESolver solverAmerican(gAmericanLeftFunc, gEuropeanRight, fOption, 0.0, S0, K, T, q, r,
                                                 sigma,
                                                 pow(M, i), alphatemp);
                solverAmerican.setUp();
                std::cout << M << "\t" << alphatemp << "\t" << solverAmerican.getAlpha() << "\t "
                          << solverAmerican.getN() << "\t" <<
                          solverAmerican.get_xLeft() << "\t" << solverAmerican.get_xRight() << "\t"
                          << solverAmerican.getXCompute() << "\t" << solverAmerican.get_tFinal() <<
                          "\t" << solverAmerican.getDeltaTau() << "\t" << solverAmerican.getDeltaX() << std::endl;
            }
        }

        void Question2Part2() {

            fOption fOption(sigma, S0, q, K, T, r);
            gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
            gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
            gAmericanRight gRightOption(sigma, S0, q, K, T, r);
            gEuropeanRight gEuropeanRight(sigma, S0, q, K, T, r);

            AmericanPDESolver solverAmerican(gAmericanLeftFunc, gEuropeanRight, fOption, 0.0, S0, K, T, q, r,
                                             sigma,
                                             M, alphatemp);
            solverAmerican.setUp();
            MatrixXd fEulerResult = solverAmerican.forwardEuler();
            std::cout << fEulerResult << std::endl;
            //std::cout << gAmericanLeftFunc.b << std::endl;
        }

        void Question2Part3() {
            fOption fOption(sigma, S0, q, K, T, r);
            gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
            gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
            gAmericanRight gAmericanRight(sigma, S0, q, K, T, r);
            // gEuropeanRight gEuropeanRight(sigma, S0, q, K, T, r);

            AmericanPDESolver solverAmerican(gAmericanLeftFunc, gAmericanRight, fOption, 0.0, S0, K, T, q, r,
                                             sigma,
                                             M, alphatemp);
            solverAmerican.setUp();
            MatrixXd fEulerResult = solverAmerican.forwardEuler();

            //gamma is printing V(i's,0)
            //solver.VApprox1 is printing V(i's,deltas)

            double gamma = solverAmerican.calculateGamma(fEulerResult); // needed to print V(i's)
            double delta = solverAmerican.calculateDelta(fEulerResult);
            double theta = solverAmerican.calculateTheta(fEulerResult);
            std::cout << "Vapprox is " << solverAmerican.calculateVApprox1(fEulerResult) << std::endl;
            std::cout << "delta is:" << delta << std::endl;
            std::cout << "gamma is:" << gamma << std::endl;
            std::cout << "theta is:" << theta << std::endl;

        }

        void Question2Part4() {
            fOption fOption(sigma, S0, q, K, T, r);
            gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
            gEuropeanRight gRightEuropeanFunc(sigma, S0, q, K, T, r);

            EuropeanPDESolver solverEuropean(gLeftEuropeanFunc, gRightEuropeanFunc, fOption, 0.0, S0, K, T, q, r,
                                             sigma,
                                             M, alphatemp);
            solverEuropean.setUp();
            MatrixXd fEulerResult = solverEuropean.forwardEuler();

            //gamma is printing V(i's,0)
            //solver.VApprox1 is printing V(i's,deltas)

            double gamma = solverEuropean.calculateGamma(fEulerResult); // needed to print V(i's)
            double delta = solverEuropean.calculateDelta(fEulerResult);
            double theta = solverEuropean.calculateTheta(fEulerResult);
            std::cout << "Vapprox is " << solverEuropean.calculateVApprox1(fEulerResult) << std::endl;
            std::cout << "delta is:" << delta << std::endl;
            std::cout << "gamma is:" << gamma << std::endl;
            std::cout << "theta is:" << theta << std::endl;
        }

        void Question2Part5() {
            BlackScholesPutOption blackScholesPutOption(S0, K, T, q, r, sigma);
            double BS_price = blackScholesPutOption.price();
            double BS_delta = blackScholesPutOption.delta();

            fOption fOption(sigma, S0, q, K, T, r);
            gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
            gAmericanRight gAmericanRight(sigma, S0, q, K, T, r);

            AmericanPDESolver solverAmerican(gAmericanLeftFunc, gAmericanRight, fOption, 0.0, S0, K, T, q, r,
                                             sigma,
                                             M, alphatemp);
            solverAmerican.setUp();
            MatrixXd fEulerResultAmerican = solverAmerican.forwardEuler();

            gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
            gEuropeanRight gRightEuropeanFunc(sigma, S0, q, K, T, r);

            EuropeanPDESolver solverEuropean(gLeftEuropeanFunc, gRightEuropeanFunc, fOption, 0.0, S0, K, T, q, r,
                                             sigma,
                                             M, alphatemp);
            solverEuropean.setUp();
            MatrixXd fEulerResultEuropean = solverEuropean.forwardEuler();
            double europeadFDPrice = solverEuropean.calculateVApprox1(fEulerResultEuropean);
            double varReductionPrice = solverAmerican.priceVarianceReduction(fEulerResultAmerican, europeadFDPrice,
                                                                             BS_price);

            std::cout << varReductionPrice << std::endl;

            //print greeks here
            std::cout << solverAmerican.calculateDelta(fEulerResultAmerican) -
                         solverEuropean.calculateDelta(fEulerResultEuropean)
                         + BS_delta << std::endl;

            //print other greeks

        }
    }
    namespace Question3 {
        long M = 1;
        double S0 = 52.0;
        double K = 55.0;
        double Lowerarrier = 45.0;
        double UpperBarrier = 62.0;
        double T = 5.0 / 12.0; // make sure to do double division
        double q = 0.005;
        double r = 0.02;
        double sigma = 0.25;
        double alphatemp = 0.5;

        void Question3part1() {
            hw8fBarrierOption fBarrierOption(S0, K, T, r, sigma, q, Lowerarrier);
            hw8gLeftBarrierOption gLeftBarrierOption(S0, K, T, r, sigma, q, Lowerarrier);
            hw8gRightBarrierOption gRightBarrierOption(S0, K, T, r, sigma, q, Lowerarrier);

            BarrierOption barrierOption(S0, K, T, q, r, sigma, Lowerarrier);

            double exactPrice = barrierOption.Price();
            MatrixXd outputMatrix = MatrixXd::Zero(4, 6);
            double AlphaTemp1 = 0.5;

            for (int i = 0; i < 4; ++i) {
                M *= 4;
                DoubleBarrierPDESolver solverAmerican(gLeftBarrierOption, gRightBarrierOption, fBarrierOption, 0.0, S0,
                                                      K, T, q, r,
                                                      sigma, AlphaTemp1, M, Lowerarrier, UpperBarrier);
                solverAmerican.setUp();
                std::cout << M << "\t" << alphatemp << "\t" << solverAmerican.getAlpha() << "\t "
                          << solverAmerican.getN() << "\t" <<
                          solverAmerican.get_xLeft() << "\t" << solverAmerican.get_xRight() << "\t"
                          << solverAmerican.getXCompute() << "\t" << solverAmerican.get_tFinal() <<
                          "\t" << solverAmerican.getDeltaTau() << "\t" << solverAmerican.getDeltaX() << std::endl;

            }

        }

    }
}





#endif //CPPCODETEST_EXAMS_HPP
