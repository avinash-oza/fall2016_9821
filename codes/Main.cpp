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


using namespace Eigen;
using namespace std;

string FILE_ROOT = "/home/avi/";
//defined at end
void writeCSVMatrix(MatrixXd &matrixToWrite, string fileName);
void hw8();
void hw9();

int main() {
    /// Keep this line to make the decimals always print out
    std::cout << std::fixed << std::setprecision(9);
	/*
	std::ofstream myfile("output1.csv");
	myfile<<(U).format(CSVFormat)<<std::endl;
	myfile.close();
	*/
//    hw8();
    hw9();
//    Question3();
//    decompositionExamples();
//    verifyCholeskyDecomposition();
//    exam2013();
//    for (int i = 10; i <= 1280; i *= 2 )
//    {
//        calculateTrinomialTreesForN(i);
//    }

    return 0;
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
//    double alphatemp = 0.45;
    double alphatemp = 5.0;

    double P_amer_bin = 4.083817051176386;

    hw8fOption fOption(sigma, S0, q, K, T, r);
    hw8gLeftOption gLeftOption(sigma, S0, q, K, T, r);
    gAmericanLeftFunc gLeftAmerican(sigma, S0, q, K, T, r);
    hw8gRightOption gRightOption(sigma, S0, q, K, T, r);

    BlackScholesOption option(S0, K, T, q, r, sigma);
    double Vexact = option.putPrice();

    for(int i = 0; i < 4; ++i) {
        M *= 4;

        EuropeanPutPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
        AmericanPutPDESolver solverAmerican(gLeftAmerican, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);

        MatrixXd fEulerEuropean = solver.CrankNicolson(SOR, tol, omega);
        double VEurApprox = solver.calculateVapprox(S0, fEulerEuropean);
//        std::cout << solverAmerican.calculateErrorPointwise(fEulerEuropean, Vexact) << std::endl;


        MatrixXd fEulerAmerican = solverAmerican.CrankNicolson(SOR, tol, omega);
//        writeCSVMatrix(fEulerAmerican, "forward_euler_american_crank.csv");
//        MatrixXd fEulerAmerican = solverAmerican.forwardEuler();
//        std::cout << fEulerAmerican << std::endl;
        std::cout << M <<","
                  << solverAmerican.calculateErrorPointwise(fEulerAmerican, P_amer_bin)
                  << ","
                  << solverAmerican.calculateErrorPointwise2(fEulerAmerican, P_amer_bin)
//                  << std::endl;
//                std::cout << M <<"," << solverAmerican.calculateDelta(fEulerAmerican) << ","
                          << solverAmerican.calculateGamma(fEulerAmerican)
                          << ","
                          << solverAmerican.calculateTheta(fEulerAmerican)
//                          << std::endl;
//            std::cout << M
                      << ","
//                      << solverAmerican.calculateVapprox(S0, fEulerAmerican)
//                      << ","
//                      << VEurApprox
//                      << ","
//                      << Vexact
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

    BlackScholesOption option(S0, K, T, q, r, sigma);
    double Vexact = option.putPrice();

    for(int i = 0; i < 4; ++i) {
        M *= 4;

        EuropeanPutPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
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


