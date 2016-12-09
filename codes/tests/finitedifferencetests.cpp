//
// Created by avi on 12/8/16.
//
#include <iostream>

#include "../EuropeanPDESolver.hpp"
#include "gtest/gtest.h"

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

    hw8fOption fOption(sigma, S0, q, K, T, r);
    hw8gLeftOption gLeftOption(sigma, S0, q, K, T, r);
    hw8gRightOption gRightOption(sigma, S0, q, K, T, r);

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
