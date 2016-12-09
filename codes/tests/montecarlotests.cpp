//
// Created by avi on 12/9/16.
//
#include <iostream>

#include "gtest/gtest.h"
#include "../RandomNumberGenerator.hpp"
#include "../MonteCarlo.hpp"
#include "../BlackScholes.hpp"

double TOL3 = std::pow(10, -10);

TEST(MonteCarloTests, MonteCarloTests_AntitheticVariables_Test)
{

    // Calculating the option price.
    double spot, strike, interest, vol, maturity, div;
    spot = 50; strike = 55; interest = 0.04; vol = 0.3; maturity = 0.5; div = 0;
    BlackScholesPutOption option(spot, strike, maturity, div, interest, vol);
    double price;

    price = option.price();

    // Creating the vector containing the value for the number of random variables.
//    int N_vector_size = 10;
    int N_vector_size = 1;
    Eigen::VectorXi N_vector = Eigen::VectorXi::Zero(N_vector_size);
//    N_vector << 1, 2, 4, 8, 16, 32, 64, 128, 256, 512;
    N_vector << 256;
    N_vector *= 10000;


    // Part a: Iverse Transform Method
    BoxMullerMethod inverseTransformMethod;
    LinearCongruentialGenerator uniformMethod;
    AntitheticMonteCarloMethod monteCarloPricer;

    MatrixXd results = monteCarloPricer.runMonteCarloForPaths(spot, strike, interest, vol, div, maturity, N_vector, inverseTransformMethod, uniformMethod, price);
    ASSERT_NEAR(results(0,2), 6.61710048461109, TOL3);


}


TEST(MonteCarloTests, MonteCarloTests_MomentMatchingAndControlVariates_Test)
{

    // Calculating the option price.
    double spot, strike, interest, vol, maturity, div;
    spot = 50; strike = 55; interest = 0.04; vol = 0.3; maturity = 0.5; div = 0;
    BlackScholesPutOption option(spot, strike, maturity, div, interest, vol);
    double price;

    price = option.price();

    // Creating the vector containing the value for the number of random variables.
//    int N_vector_size = 10;
    int N_vector_size = 1;
    Eigen::VectorXi N_vector = Eigen::VectorXi::Zero(N_vector_size);
//    N_vector << 1, 2, 4, 8, 16, 32, 64, 128, 256, 512;
    N_vector << 256;
    N_vector *= 10000;


    // Part a: Iverse Transform Method
    BoxMullerMethod inverseTransformMethod;
    LinearCongruentialGenerator uniformMethod;
    MomentMatchingAndControlVariateMonteCarloMethod monteCarloPricer;

    MatrixXd results = monteCarloPricer.runMonteCarloForPaths(spot, strike, interest, vol, div, maturity, N_vector, inverseTransformMethod, uniformMethod, price);
    ASSERT_NEAR(results(0,2), 6.61812397917852, TOL3);




}
