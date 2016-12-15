//
// Created by Tushar Chawla on 12/9/16.
//

#ifndef CPPCODETEST_EXAMS_HPP
#define CPPCODETEST_EXAMS_HPP

#include <iostream>
#include "BinomialTrees.hpp"
#include "AmericanPDESolver.hpp"
#include "BarrierOption.hpp"
#include "BarrierPDE.hpp"
#include "MonteCarlo.hpp"
#include <fstream>

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

            //this will print tree as well. prinitng for both american and european options is in european class
            //in BinomialTree.hpp file on page 86.


        }

        void Question1Part2() {


            AmericanBinomialTreePricer americanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            BlackScholesPutOption blackScholesPutOption(S, K, T, q, r, sigma);

            std::cout << americanBinomialTreePricer.calculateOptimalN(1, 1000, 0.0001) << std::endl;
        }

        void Question1Part3() {
            //check the value of optimal N from above
            long optimalN = 191;
            EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            TREE_RESULT european_tree_result = europeanBinomialTreePricer.calculateTree(optimalN);
            AmericanBinomialTreePricer americanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
            TREE_RESULT america_tree_result = americanBinomialTreePricer.calculateTree(optimalN);
            BlackScholesPutOption blackScholesPutOption(S, K, T, q, r, sigma);

            std::cout << europeanBinomialTreePricer.extractPrice(european_tree_result) << std::endl;
            // std::cout << europeanBinomialTreePricer.extractDelta(european_tree_result) << std:: endl;
            //std::cout << europeanBinomialTreePricer.extractGamma(european_tree_result)<< std::endl;
            //std::cout << europeanBinomialTreePricer.extractTheta(european_tree_result)<< std:: endl;


            /* std::cout << "Now printing the variance reduction American call price " << std::endl;
             std::cout << americanBinomialTreePricer.varianceReductionPrice(america_tree_result,european_tree_result,blackScholesPutOption) << std:: endl;
             std::cout << americanBinomialTreePricer.varianceReductionDelta(america_tree_result,european_tree_result,blackScholesPutOption) << std:: endl;
             std::cout << americanBinomialTreePricer.varianceReductionGamma(america_tree_result,european_tree_result,blackScholesPutOption) << std:: endl;
             std::cout << americanBinomialTreePricer.varianceReductionTheta(america_tree_result,european_tree_result,blackScholesPutOption) << std:: endl;
             */


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
            //gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
            gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
            gAmericanRight gRightOption(sigma, S0, q, K, T, r);
            //gEuropeanRight gEuropeanRight(sigma, S0, q, K, T, r);


            for (int i = 1; i <= 4; ++i) {

                AmericanPDESolver solverAmerican(gAmericanLeftFunc, gRightOption, fOption, 0.0, S0, K, T, q, r,
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
        void Question2Part3()
        {
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

            //Printing of V(i's) are in EuropeanPDESolver.hpp on line 150


            solverAmerican.printFDTree(fEulerResult);

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
            double BS_gamma = blackScholesPutOption.gamma();
            double BS_theta = blackScholesPutOption.theta();

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

            double europeadFDPrice = solverEuropean.calculateVApprox1(fEulerResultEuropean); //in here use VApprox1
            double varReductionPrice = solverAmerican.priceVarianceReduction(fEulerResultAmerican, europeadFDPrice,
                                                                             BS_price);


            std::cout << varReductionPrice << std::endl;

            std::cout << "print greeks here" << std::endl;

            std::cout << solverAmerican.calculateDelta(fEulerResultAmerican) -
                         solverEuropean.calculateDelta(fEulerResultEuropean)
                         + BS_delta << std::endl;

            std::cout << solverAmerican.calculateGamma(fEulerResultAmerican) - //turn off printing in gamma function
                         solverEuropean.calculateGamma(fEulerResultEuropean)
                         + BS_gamma << std::endl;

            std::cout << solverAmerican.calculateTheta(fEulerResultAmerican) -
                         solverEuropean.calculateTheta(fEulerResultEuropean)
                         + BS_theta << std::endl;

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

        // The functions
        hw8fBarrierOption fBarrierOption(S0, K, T, r, sigma, q, Lowerarrier);
        hw8gLeftBarrierOption gLeftBarrierOption(S0, K, T, r, sigma, q, Lowerarrier);
        hw8gRightBarrierOption gRightBarrierOption(S0, K, T, r, sigma, q, Lowerarrier);

        BarrierOption barrierOption(S0, K, T, q, r, sigma, Lowerarrier);

        void Question3Part1() {

            double exactPrice = barrierOption.Price();
            MatrixXd outputMatrix = MatrixXd::Zero(4, 10);
            double AlphaTemp1 = 0.5;

            for (int i = 0; i < 4; ++i) {
                M *= 4;
                DoubleBarrierPDESolver solverAmerican(gLeftBarrierOption, gRightBarrierOption, fBarrierOption, 0.0, S0,
                                                      K, T, q, r,
                                                      sigma, AlphaTemp1, M, Lowerarrier, UpperBarrier);

                solverAmerican.setUp();
                outputMatrix(i, 0) = M;
                outputMatrix(i, 1) = alphatemp;
                outputMatrix(i, 2) = solverAmerican.getAlpha();
                outputMatrix(i, 3) = solverAmerican.getN();
                outputMatrix(i, 4) = solverAmerican.get_xLeft();
                outputMatrix(i, 5) = solverAmerican.get_xRight();
                outputMatrix(i, 6) = solverAmerican.getXCompute();
                outputMatrix(i, 7) = solverAmerican.get_tFinal();
                outputMatrix(i, 8) = solverAmerican.getDeltaTau();
                outputMatrix(i, 9) = solverAmerican.getDeltaX();
            }
            std::cout << outputMatrix << std::endl;
        }
        // For this part reset M = 4

        void Question3Part2() {
            M = 4;
            double AlphaTemp2 = 0.5;
            DoubleBarrierPDESolver solverAmerican(gLeftBarrierOption, gRightBarrierOption, fBarrierOption, 0.0, S0,
                                                  K, T, q, r,
                                                  sigma, AlphaTemp2, M, Lowerarrier, UpperBarrier);
            solverAmerican.setUp();

            MatrixXd outputMatrix = solverAmerican.forwardEuler();
            std::cout << outputMatrix << std::endl;
        }

        void Question3Part3() {
            M = 4;
            double AlphaTemp3 = 0.5;
            DoubleBarrierPDESolver solverAmerican(gLeftBarrierOption, gRightBarrierOption, fBarrierOption, 0.0, S0,
                                                  K, T, q, r,
                                                  sigma, AlphaTemp3, M, Lowerarrier, UpperBarrier);
            solverAmerican.setUp();

            MatrixXd forwardSolution = solverAmerican.forwardEuler();

            double Price, Delta, Gamma, Theta;
            Price = solverAmerican.calculateVApprox(forwardSolution);
            std::cout << Price << std::endl;

            Delta = solverAmerican.calculateDelta(forwardSolution);
            std::cout << Delta << std::endl;

            Gamma = solverAmerican.calculateGamma(forwardSolution);
            std::cout << Gamma << std::endl;

            Theta = solverAmerican.calculateTheta(forwardSolution);
            std::cout << Theta << std::endl;
        }
    }

}






#endif //CPPCODETEST_EXAMS_HPP
