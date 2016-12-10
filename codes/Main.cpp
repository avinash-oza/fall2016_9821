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
#include "RandomNumberGenerator.hpp"
#include "MonteCarlo.hpp"
//#include "BarrierPDE.hpp"
#include "Exams.hpp"
#include "BarrierPDE.hpp"

using namespace Eigen;
using namespace std;

string FILE_ROOT = "/home/avi/";
//defined at end
void writeCSVMatrix(MatrixXd &matrixToWrite, string fileName);
void hw8();
void hw9();
void Question3();
void hw9Barrier();
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
//    hw10();
    hw9Barrier();
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


void hw9Barrier()
{
    // Part 1 : Domain Discretization  M = {4, 16, 64, 256} alpha = 0.4
    //std::cout << "Domain Discretization: alpha_temp = 0.4\n\n";
    double S0 = 42;
    double K = 40;
    double q = 0.03;
    double Sigma = 0.3;
    double B = 35;
    double r = 0.05;
    double T = 0.5;

    double tol = std::pow(10, -6);
    double Omega = 1.2;
    long M = 1;

    hw8fBarrierOption fBarrierOption(S0, K, T, r, Sigma, q, B);
    hw8gLeftBarrierOption gLeftBarrierOption(S0, K, T, r, Sigma, q, B);
    hw8gRightBarrierOption gRightBarrierOption(S0, K, T, r, Sigma, q, B);

    BarrierOption barrierOption(S0, K, T, q, r, Sigma, B);
    double exactPrice = barrierOption.Price();
    MatrixXd outputMatrix = MatrixXd::Zero(4, 6);
    double AlphaTemp1 = 0.4;

    for (int i = 0; i < 4; ++i)
    {
        M *= 4;
    	DownOutCallPDESolver Barrier1(gLeftBarrierOption, gRightBarrierOption, fBarrierOption, 0.0, S0, K, T, q, r,
                                      Sigma, AlphaTemp1, M, B);
    	Barrier1.setUp();
        long N_left = Barrier1.get_N_left();
        std::cout << M << std::endl;

        auto ApproxMatrix = Barrier1.CrankNicolson(SOR, tol, Omega);
        outputMatrix(i, 0) = ApproxMatrix(M, N_left);	// The value of the option is in last raw of each, first column of each approximation matrix
        outputMatrix(i, 1) = Barrier1.calculateVApprox1(ApproxMatrix);
        outputMatrix(i, 2) = Barrier1.calculateErrorPointwise(ApproxMatrix, exactPrice);
        outputMatrix(i, 3) = Barrier1.calculateDelta(ApproxMatrix);
        outputMatrix(i, 4) = Barrier1.calculateGamma(ApproxMatrix);
        outputMatrix(i, 5) = Barrier1.calculateTheta(ApproxMatrix);
    }
    std:cout << outputMatrix << std::endl;
}

void hw10()
{
    double S = 41.0;
    double K = 40.0;
    double T = 1.0; // make sure to do double division
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;
    OPTION_TYPE optionType = PUT;
    double B = 45.0;
    BARRIER_TYPE barrierType = DOWN_AND_OUT;
    //MAKE SURE TO CHANGE PAYOFF FUNCTION
    AmericanBinomialTreePricer americanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
    EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);
    BlackScholesPutOption blackScholesPutOption(S, K, T, q, r, sigma);

    TREE_RESULT pricerResult = americanBinomialTreePricer.calculateTree(20);
    TREE_RESULT europeanPricerResult = europeanBinomialTreePricer.calculateTree(20);

//    TREE_RESULT pricerResult2 = barrierTrinomialTreePricer.calculateTree(1000);

//    std::cout << blackScholesPutOption.gamma() << std::endl;

    std::cout << europeanBinomialTreePricer.varianceReductionGamma(pricerResult, europeanPricerResult, blackScholesPutOption)
//            <<","
//                << barrierTrinomialTreePricer.extractPrice(pricerResult2)
//               <<","
//                << blackScholesPutOption.Price()
              << std::endl;
}


void hw9()
{
    long M = 64;
    double tol = std::pow(10, -6);
    double omega = 1.2;

    double sigma = 0.35;
    double S0 = 41;
    double q = 0.02;
    double K = 40;
    double T = 0.75;
    double r = 0.04;
    double alphatemp = 0.45;

    double P_amer_bin = 4.083817051176386;

    fOption fOption(sigma, S0, q, K, T, r);
    gEuropeanLeft gLeftEuropeanFunc(sigma, S0, q, K, T, r);
    gAmericanLeftFunc gAmericanLeftFunc(sigma, S0, q, K, T, r);
    gAmericanRight gRightOption(sigma, S0, q, K, T, r);
    gEuropeanRight gEuropeanRight(sigma, S0, q, K, T, r);

    BlackScholesPutOption option(S0, K, T, q, r, sigma);
    double Vexact = option.price();

    EuropeanPDESolver solver(gLeftEuropeanFunc, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
    solver.setUp();
    AmericanPDESolver solverAmerican(gAmericanLeftFunc, gEuropeanRight, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
    solverAmerican.setUp();

    MatrixXd fEulerResultEuropean = solver.forwardEuler();
    double fEulerEuropeanResult = solver.calculateVapprox(fEulerResultEuropean);

    MatrixXd fEulerResult = solverAmerican.forwardEuler();
    MatrixXd CNResult = solverAmerican.CrankNicolson(SOR, tol, omega);

    double forwardEulerVApproxEuropean = solver.calculateVapprox(fEulerResultEuropean);
            std::cout << M
                      << ","
                      << solverAmerican.calculateVapprox(fEulerResult)
                      << ","
//                      << VEurApprox
//                      << ","
                      << Vexact
                      << solverAmerican.priceVarianceReduction(fEulerResult, forwardEulerVApproxEuropean, Vexact)
                      << ","
                      << solverAmerican.calculateErrorPointWiseVarianceReduction(fEulerResult, forwardEulerVApproxEuropean, Vexact, P_amer_bin)
                      << std::endl;

//    writeCSVMatrix(fEulerAmerican, "forward_euler_american.csv");

//    writeCSVMatrix(fEulerResult, "/home/avi/forwardEuler.csv");
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

    fOption fOption(sigma, S0, q, K, T, r);
    gEuropeanLeft gLeftOption(sigma, S0, q, K, T, r);
    gAmericanRight gRightOption(sigma, S0, q, K, T, r);

    BlackScholesPutOption option(S0, K, T, q, r, sigma);
    double Vexact = option.price();

    for(int i = 0; i < 4; ++i) {
        M *= 4;

        EuropeanPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
        solver.setUp();
        MatrixXd fEulerResult = solver.CrankNicolson(SOR, tol, omega);
    std::cout << M << std::endl;
        std::cout << solver.calculateErrorPointwise(fEulerResult, Vexact) << std::endl;
//        std::cout << solver.calculateErrorPointwise2(fEulerResult, Vexact) << std::endl;
    }
//    writeCSVMatrix(fEulerResult, "/home/avi/forwardEuler.csv");
}

void Question3()
{
    // Calculating the option price.
    double spot, strike, interest, vol, maturity, div;
    spot = 50; strike = 55; interest = 0.04; vol = 0.3; maturity = 0.5; div = 0;
    BlackScholesPutOption option(spot, strike, maturity, div, interest, vol);
    double price;

    price = option.price();

    // Creating the vector containing the value for the number of random variables.
    int N_vector_size = 10;
//    int N_vector_size = 1;
    Eigen::VectorXi N_vector = Eigen::VectorXi::Zero(N_vector_size);
    N_vector << 1, 2, 4, 8, 16, 32, 64, 128, 256, 512;
//    N_vector << 1;
    N_vector *= 10000;


    // Part a: Iverse Transform Method
    BoxMullerMethod inverseTransformMethod;
    LinearCongruentialGenerator uniformMethod;
    AntitheticMonteCarloMethod monteCarloPricer;

    MatrixXd results = monteCarloPricer.runMonteCarloForPaths(spot, strike, interest, vol, div, maturity, N_vector, inverseTransformMethod, uniformMethod, price);
    std::cout << results << std::endl;
}

void writeCSVMatrix(MatrixXd &matrixToWrite, string fileName)
{
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ",", "\n");
    std::ofstream myfile(FILE_ROOT + fileName);
    myfile<<(matrixToWrite).format(CSVFormat) << std::endl;
    myfile.close();
    std::cout << "Wrote output to " << fileName << std::endl;
}


