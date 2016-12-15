//
// Created by avi on 12/8/16.
//

#include <iostream>

#include "../BinomialTrees.hpp"
#include "../TrinomialTrees.hpp"
#include "gtest/gtest.h"

double TOL = std::pow(10, -6);

TEST(BinomialTest, BinomialTest_BarrierCallTest_Test)
{
    double S = 50.0;
    double K = 48.0;
    double T = 8.0/12.0; // make sure to do double division
    double q = 0.01;
    double r = 0.02;
    double sigma = 0.3;
    double B = 45.0;
    OPTION_TYPE optionType = CALL;
    BARRIER_TYPE barrierType = DOWN_AND_OUT;

    BarrierOptionBinomialTreePricer barrierOptionBinomialTreePricer(S, K, T, q, r, sigma, B, barrierType, optionType);

    TREE_RESULT oneResult  = barrierOptionBinomialTreePricer.calculateTree(1000); // calculate for 1 iteration
    double result = barrierOptionBinomialTreePricer.extractPrice(oneResult);
    EXPECT_NEAR(result , 4.295672979888, TOL );

}


TEST(BinomialTest, BinomialTest_BinomialTest_AmericanOption_Test)
{

    double S = 41.0;
    double K = 40.0;
    double T = 1.0; // make sure to do double division
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;
    OPTION_TYPE optionType = PUT;

    AmericanBinomialTreePricer americanBinomialTreePricer(S, K, T, q, r, sigma, optionType);

    TREE_RESULT oneResult  = americanBinomialTreePricer.calculateTree(20); // calculate for 1 iteration

    EXPECT_NEAR(americanBinomialTreePricer.extractPrice(oneResult) , 3.979843704, TOL );
    EXPECT_NEAR(americanBinomialTreePricer.extractDelta(oneResult) , -0.3895857885, TOL );
    EXPECT_NEAR(americanBinomialTreePricer.extractGamma(oneResult) , 0.03283395365, TOL );
    EXPECT_NEAR(americanBinomialTreePricer.extractTheta(oneResult) , -2.050180819, TOL );


    oneResult  = americanBinomialTreePricer.BlackScholesWithRichardsonExtrapolation(20); // calculate for 1 iteration

    EXPECT_NEAR(americanBinomialTreePricer.extractPrice(oneResult) , 3.976027001762, TOL );
    EXPECT_NEAR(americanBinomialTreePricer.extractDelta(oneResult) , -0.388997632963443, TOL );
    EXPECT_NEAR(americanBinomialTreePricer.extractGamma(oneResult) , 0.032005666610067, TOL );
    EXPECT_NEAR(americanBinomialTreePricer.extractTheta(oneResult) , -1.9825255477319, TOL );

}



TEST(BinomialTest, BinomialTest_EuropeanOption_Test)
{

    double S = 41.0;
    double K = 40.0;
    double T = 1.0; // make sure to do double division
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;
    OPTION_TYPE optionType = PUT;

    EuropeanBinomialTreePricer europeanBinomialTreePricer(S, K, T, q, r, sigma, optionType);

    TREE_RESULT oneResult  = europeanBinomialTreePricer.calculateTree(20); // calculate for 1 iteration

    EXPECT_NEAR(europeanBinomialTreePricer.extractPrice(oneResult) , 3.9092991120000, TOL );
    EXPECT_NEAR(europeanBinomialTreePricer.extractDelta(oneResult) , -0.3797719590000, TOL );
    EXPECT_NEAR(europeanBinomialTreePricer.extractGamma(oneResult) , 0.0314586080000, TOL );
    EXPECT_NEAR(europeanBinomialTreePricer.extractTheta(oneResult) , -1.9560716110000, TOL );


    oneResult  = europeanBinomialTreePricer.BlackScholesWithRichardsonExtrapolation(20); // calculate for 1 iteration

    EXPECT_NEAR(europeanBinomialTreePricer.extractPrice(oneResult) , 3.900063684000, TOL );
    EXPECT_NEAR(europeanBinomialTreePricer.extractDelta(oneResult) , -0.378699586000, TOL );
    EXPECT_NEAR(europeanBinomialTreePricer.extractGamma(oneResult) , 0.030610531000, TOL );
    EXPECT_NEAR(europeanBinomialTreePricer.extractTheta(oneResult) , -1.887709740000, TOL );
    EXPECT_EQ(europeanBinomialTreePricer.calculateOptimalN(1, 1000, 0.0001), 331);

}


////////////////////////////// TRINOMIAL TESTS/////////////////////////////


TEST(TrinomialTest, TrinomialTest_AmericanTrinomial_Test)
{
    double S = 41.0;
    double K = 40.0;
    double T = 1.0; // make sure to do double division
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.3;
    OPTION_TYPE optionType = PUT;

    AmericanTrinomialTreePricer americanTrinomialTreePricer(S, K, T, q, r, sigma, optionType);

    TREE_RESULT oneResult  = americanTrinomialTreePricer.calculateTree(40); // calculate for 1 iteration


    EXPECT_NEAR(americanTrinomialTreePricer.extractPrice(oneResult) , 3.972490588, TOL );
    EXPECT_NEAR(americanTrinomialTreePricer.extractDelta(oneResult) , -0.387651998, TOL );
    EXPECT_NEAR(americanTrinomialTreePricer.extractGamma(oneResult) , 0.03204947773, TOL );
    EXPECT_NEAR(americanTrinomialTreePricer.extractTheta(oneResult) , -2.00105521, TOL );

    oneResult  = americanTrinomialTreePricer.BlackScholesWithRichardsonExtrapolation(40); // calculate for 1 iteration

    EXPECT_NEAR(americanTrinomialTreePricer.extractPrice(oneResult) , 3.970801179, TOL );
    EXPECT_NEAR(americanTrinomialTreePricer.extractDelta(oneResult) , -0.388665148, TOL );
    EXPECT_NEAR(americanTrinomialTreePricer.extractGamma(oneResult) , 0.032085761, TOL );
    EXPECT_NEAR(americanTrinomialTreePricer.extractTheta(oneResult) , -1.988626278, TOL );

}