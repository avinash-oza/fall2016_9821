//
// Created by avi on 12/8/16.
//
#include <iostream>

#include "../EuropeanPDESolver.hpp"
#include "gtest/gtest.h"
#include "../AmericanPDESolver.hpp"
#include "../BarrierOption.hpp"
#include "../BarrierPDE.hpp"

double TOL2 = std::pow(10, -9);

TEST(FiniteDifferenceTests, FiniteDifferenceTests_EuropeanPutPricing_Test)
{

    hw8gLeft  gLeft;
    hw8gRight gRight;
    hw8f f;
    uExact uExact1;

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

    fOption fOption(sigma, S0, q, K, T, r);
    gEuropeanLeft gLeftOption(sigma, S0, q, K, T, r);
    gAmericanRight gRightOption(sigma, S0, q, K, T, r);

    BlackScholesPutOption option(S0, K, T, q, r, sigma);
    double Vexact = option.price();

    EuropeanPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
    solver.setUp();
    MatrixXd fEulerResult = solver.forwardEuler();
    MatrixXd bEulerResult = solver.backwardEuler(LU, tol, omega);
    MatrixXd CNResult = solver.CrankNicolson(SOR, tol, omega);

    // forward euler values
    ASSERT_NEAR(solver.calculateErrorPointwise(fEulerResult, Vexact), 0.016102244, TOL2);
    ASSERT_NEAR(solver.calculateErrorPointwise2(fEulerResult, Vexact), 0.015378482, TOL2);
    ASSERT_NEAR(solver.calculateDelta(fEulerResult), -0.371230484, TOL2);
    ASSERT_NEAR(solver.calculateGamma(fEulerResult), 0.029584544, TOL2);
    ASSERT_NEAR(solver.calculateTheta(fEulerResult), -2.65551477, TOL2);

    // forward euler values
    ASSERT_NEAR(solver.calculateErrorPointwise(bEulerResult, Vexact), 0.001709905, TOL2);
    ASSERT_NEAR(solver.calculateDelta(bEulerResult), -0.371022935, TOL2);
    ASSERT_NEAR(solver.calculateGamma(bEulerResult), 0.02992175, TOL2);
    ASSERT_NEAR(solver.calculateTheta(bEulerResult), -2.668583705, TOL2);


    // crank nicolson values
    ASSERT_NEAR(solver.calculateErrorPointwise(CNResult, Vexact), 0.007221706, TOL2);
    ASSERT_NEAR(solver.calculateErrorPointwise2(CNResult, Vexact), 0.006498, TOL2);
    ASSERT_NEAR(solver.calculateDelta(CNResult), -0.371129896, TOL2);
    ASSERT_NEAR(solver.calculateGamma(CNResult), 0.029750817, TOL2);
    ASSERT_NEAR(solver.calculateTheta(CNResult), -2.661998976, TOL2);

}


TEST(FiniteDifferenceTests, FiniteDifferenceTests_AmericanPutPricing_Test)
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

        MatrixXd fEulerEuropean = solver.forwardEuler();
        MatrixXd fEulerAmerican = solverAmerican.forwardEuler();

        MatrixXd CNResultAmerican = solverAmerican.CrankNicolson(SOR, tol, omega);
        MatrixXd CNResultEuropean = solver.CrankNicolson(SOR, tol, omega);

        double forwardEulerVApproxEuropean = solver.calculateVApprox(fEulerEuropean);
        double CNVApproxEuropean = solver.calculateVApprox(CNResultEuropean);

    // forward euler values
    ASSERT_NEAR(solverAmerican.calculateErrorPointwise(fEulerAmerican, P_amer_bin), 0.017420108129288, TOL2);
    ASSERT_NEAR(solverAmerican.calculateErrorPointwise2(fEulerAmerican, P_amer_bin), 0.016680703725282, TOL2);
    ASSERT_NEAR(solverAmerican.calculateDelta(fEulerAmerican), -0.379022202454181, TOL2);
    ASSERT_NEAR(solverAmerican.calculateGamma(fEulerAmerican), 0.030667406086139, TOL2);
    ASSERT_NEAR(solverAmerican.calculateTheta(fEulerAmerican), -2.76313760830757, TOL2);
    ASSERT_NEAR(solverAmerican.priceVarianceReduction(fEulerAmerican, forwardEulerVApproxEuropean, Vexact), 4.08513491516998, TOL2);
    ASSERT_NEAR(solverAmerican.calculateErrorPointWiseVarianceReduction(fEulerAmerican, forwardEulerVApproxEuropean, Vexact, P_amer_bin), 0.00131786399359, TOL2);


    // crank nicolson values
    ASSERT_NEAR(solverAmerican.calculateErrorPointwise(CNResultAmerican, P_amer_bin), 0.005435222342801, TOL2);
    ASSERT_NEAR(solverAmerican.calculateErrorPointwise2(CNResultAmerican, P_amer_bin), 0.004696240280308, TOL2);
    ASSERT_NEAR(solverAmerican.calculateDelta(CNResultAmerican), -0.378724733614418, TOL2);
    ASSERT_NEAR(solverAmerican.calculateGamma(CNResultAmerican), 0.030833522375094, TOL2);
    ASSERT_NEAR(solverAmerican.calculateTheta(CNResultAmerican), -2.77039243739341, TOL2);
    ASSERT_NEAR(solverAmerican.priceVarianceReduction(CNResultAmerican, CNVApproxEuropean, Vexact), 4.08203056760127, TOL2);
    ASSERT_NEAR(solverAmerican.calculateErrorPointWiseVarianceReduction(CNResultAmerican, CNVApproxEuropean, Vexact, P_amer_bin), 0.001786483575114, TOL2);

}


TEST(FiniteDifferenceTests, FiniteDifferenceTests_DownAndOutBarrierOption_Test)
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
    long M = 16;

    hw8fBarrierOption fBarrierOption(S0, K, T, r, Sigma, q, B);
    hw8gLeftBarrierOption gLeftBarrierOption(S0, K, T, r, Sigma, q, B);
    hw8gRightBarrierOption gRightBarrierOption(S0, K, T, r, Sigma, q, B);

    BarrierOption barrierOption(S0, K, T, q, r, Sigma, B);
    double exactPrice = barrierOption.Price();
    MatrixXd outputMatrix = MatrixXd::Zero(4, 6);
    double AlphaTemp1 = 0.4;

    DownOutCallPDESolver downOutCallPDESolver(gLeftBarrierOption, gRightBarrierOption, fBarrierOption, 0.0, S0, K, T, q, r,
                                  Sigma, AlphaTemp1, M, B);
    downOutCallPDESolver.setUp();
    long N_left = downOutCallPDESolver.get_N_left();
    std::cout << M << std::endl;

    auto fEulerAmerican = downOutCallPDESolver.forwardEuler();
    auto fEulerCN = downOutCallPDESolver.CrankNicolson(SOR, tol, Omega);

    // forward euler values
    ASSERT_NEAR(downOutCallPDESolver.calculateDelta(fEulerAmerican), 0.707288582514209, TOL2);
    ASSERT_NEAR(downOutCallPDESolver.calculateGamma(fEulerAmerican), 0.024436109375867, TOL2);
    ASSERT_NEAR(downOutCallPDESolver.calculateTheta(fEulerAmerican), -2.45879570557505, TOL2);
    ASSERT_NEAR(downOutCallPDESolver.calculateVApprox1(fEulerAmerican), 4.46444038054815, TOL2);


    // crank nicolson values
    ASSERT_NEAR(downOutCallPDESolver.calculateDelta(fEulerCN), 0.707582938228157, TOL2);
    ASSERT_NEAR(downOutCallPDESolver.calculateGamma(fEulerCN), 0.025856537817444, TOL2);
    ASSERT_NEAR(downOutCallPDESolver.calculateTheta(fEulerCN), -2.50666267382413, TOL2);
    ASSERT_NEAR(downOutCallPDESolver.calculateVApprox1(fEulerCN), 4.42851597070506, TOL2);

}





