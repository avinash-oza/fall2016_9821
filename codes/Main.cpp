#include <tuple>
#include <iomanip>
#include<iostream>
#include<fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include "Eigen/Dense"
#include "IterativeMethods.hpp"
#include "StoppingCriterionSolvers.hpp"
#include "Solvers.hpp"
#include "OptionPricers.hpp"
#include "BlackScholes.hpp"
#include "PDESolver.hpp"
#include "EuropeanPDESolver.hpp"
#include "AmericanPDESolver.hpp"
#include "BinomialTrees.hpp"
#include "TrinomialTrees.hpp"
#include "BarrierOption.hpp"


using namespace Eigen;
using namespace std;

string FILE_ROOT = "/home/avi/";
//defined at end
void writeCSVMatrix(MatrixXd &matrixToWrite, string fileName);
void hw8();
void hw9();
void hw10();

int main() {
    /// Keep this line to make the decimals always print out
    std::cout << std::fixed << std::setprecision(12);
	/*
	std::ofstream myfile("output1.csv");
	myfile<<(U).format(CSVFormat)<<std::endl;
	myfile.close();
	*/
//    hw8();
//    hw9();
    hw10();
//    Question3();
//    decompositionExamples();
//    verifyCholeskyDecomposition();
//    exam2013();
//    double S = 41.0;
//    double K = 40.0;
//    double T = 1.0;
//    double q = 0.01;
//    double r = 0.03;
//    double sigma = 0.3;
//
//    EuropeanTrinomialTreePricer trinomialTreePricer(S, K, T, q, r, sigma);
//    for (int i = 10; i <= 1280; i *= 2 )
//    {
//
//        TREE_RESULT pricerResult = trinomialTreePricer.BlackScholesWithRichardsonExtrapolation(i);

//        std::cout << trinomialTreePricer.extractPrice(pricerResult) //- calculateTrinomialTreesForN(i)
//        std::cout << calculateTrinomialTreesForN(i)
//                << ","
//                    << trinomialTreePricer.extractDelta(pricerResult)
//                << ","
//                  << trinomialTreePricer.extractGamma(pricerResult)
//                << ","
//                  << trinomialTreePricer.extractTheta(pricerResult)
//                  << std::endl;
//    }

    return 0;
}

void hw10()
{
    double S = 50.0;
    double K = 48.0;
    double T = 8.0/12.0; // make sure to do double division
    double q = 0.01;
    double r = 0.02;
    double sigma = 0.3;
    double B = 45.0;
    //MAKE SURE TO CHANGE PAYOFF FUNCTION
    BarrierOptionBinomialTreePricer barrierOptionBinomialTreePricer(S, K, T, q, r, sigma, B, DOWN_AND_OUT, CALL);
    BarrierTrinomialTreePricer barrierTrinomialTreePricer(S, K, T, q, r, sigma, B, DOWN_AND_OUT, CALL);
    BarrierOption barrierOption(S, K, T, q, r, sigma, B);

    TREE_RESULT pricerResult = barrierOptionBinomialTreePricer.calculateTree(1000);
    TREE_RESULT pricerResult2 = barrierTrinomialTreePricer.calculateTree(1000);

    std::cout << barrierOptionBinomialTreePricer.extractPrice(pricerResult)
            <<","
                << barrierTrinomialTreePricer.extractPrice(pricerResult2)
               <<","
                << barrierOption.Price()
              << std::endl;
}


void hw9()
{
    long M = 1;
    cout << setprecision(16);
    double tol = std::pow(10, -6);
    double omega = 1.2;

    double sigma = 0.35;
    double S0 = 41;
    double q = 0.02;
    double K = 40;
    double T = 0.75;
    double r = 0.04;
    double alphatemp = 0.45;
//    double alphatemp = 5.0;

    double P_amer_bin = 4.083817051176386;

    hw8fOption fOption(sigma, S0, q, K, T, r);
    hw8gLeftOption gLeftOption(sigma, S0, q, K, T, r);
    gAmericanLeftFunc gLeftAmerican(sigma, S0, q, K, T, r);
    hw8gRightOption gRightOption(sigma, S0, q, K, T, r);

    BlackScholesPutOption option(S0, K, T, q, r, sigma);
    double Vexact = option.price();

    for(int i = 0; i < 4; ++i) {
        M *= 4;

        EuropeanPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
        solver.setUp();
        AmericanPDESolver solverAmerican(gLeftAmerican, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
        solverAmerican.setUp();

        MatrixXd fEulerEuropean = solver.CrankNicolson(SOR, tol, omega);
        double VEurApprox = solver.calculateVapprox(fEulerEuropean);
//        std::cout << solverAmerican.calculateErrorPointwise(fEulerEuropean, Vexact) << std::endl;


        MatrixXd fEulerAmerican = solverAmerican.CrankNicolson(SOR, tol, omega);
//        VectorXd t = solverAmerican.findSEarlyExerciseSoptimal(fEulerAmerican);
//        writeCSVMatrix(fEulerAmerican, "forward_euler_american_crank.csv");
//        MatrixXd fEulerAmerican = solverAmerican.forwardEuler();
//        std::cout << fEulerAmerican << std::endl;
        std::cout << M <<","
//                  << solverAmerican.calculateErrorPointwise(fEulerAmerican, P_amer_bin)
//                  << ","
//                  << solverAmerican.calculateErrorPointwise2(fEulerAmerican, P_amer_bin)
//                  << ","
//                  << std::endl;
//                std::cout << M <<","
//                    << solverAmerican.calculateDelta(fEulerAmerican)
//                    << ","
//                          << solverAmerican.calculateGamma(fEulerAmerican)
//                          << ","
//                          << solverAmerican.calculateTheta(fEulerAmerican)
                          << std::endl;
            std::cout << M
                      << ","
                      << solverAmerican.calculateVapprox(fEulerAmerican)
                      << ","
                      << VEurApprox
                      << ","
                      << Vexact
                      << solverAmerican.priceVarianceReduction(fEulerAmerican, VEurApprox, Vexact)
                      << ","
                      << solverAmerican.calculateErrorPointWiseVarianceReduction(fEulerAmerican, VEurApprox, Vexact, P_amer_bin)
                      << std::endl;

//    writeCSVMatrix(fEulerAmerican, "forward_euler_american.csv");

//    writeCSVMatrix(fEulerResult, "/home/avi/forwardEuler.csv");
    }
}

void hw8()
{
    hw8gLeft  gLeft;
    hw8gRight gRight;
    hw8f f;
    uExact uExact1;


    double xLeft = -2.0;
    double xRight = 2.0;
    double tauFinal = 1.0;
    //long M = 2;
    //long N = 8;
	long M = 1;
    cout << setprecision(16);
    double tol = std::pow(10, -6);
    double omega = 1.2;

	double sigma = 0.35;
	double S0 = 41;
	double q = 0.02;
	double K = 40;
	double T = 0.75;
	double r = 0.04;
	double alphatemp = 0.45;

    hw8fOption fOption(sigma, S0, q, K, T, r);
    hw8gLeftOption gLeftOption(sigma, S0, q, K, T, r);
    hw8gRightOption gRightOption(sigma, S0, q, K, T, r);

    BlackScholesPutOption option(S0, K, T, q, r, sigma);
    double Vexact = option.price();

    for(int i = 0; i < 4; ++i) {
        M *= 4;

        EuropeanPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
        solver.setUp();
        MatrixXd fEulerResult = solver.CrankNicolson(SOR, tol, omega);
//    std::cout << fEulerResult << std::endl;
        std::cout << solver.calculateErrorPointwise(fEulerResult, Vexact) << std::endl;
//        std::cout << solver.calculateErrorPointwise2(fEulerResult, Vexact) << std::endl;
    }
//    writeCSVMatrix(fEulerResult, "/home/avi/forwardEuler.csv");
}

void writeCSVMatrix(MatrixXd &matrixToWrite, string fileName)
{
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ",", "\n");
    std::ofstream myfile(FILE_ROOT + fileName);
    myfile<<(matrixToWrite).format(CSVFormat) << std::endl;
    myfile.close();
    std::cout << "Wrote output to " << fileName << std::endl;
}


