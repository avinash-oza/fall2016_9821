//
// Created by avi on 9/26/16.
//

#ifndef CPPCODETEST_BLACKSCHOLES_HPP
#define CPPCODETEST_BLACKSCHOLES_HPP

#include <cmath>
#include <boost/math/distributions/normal.hpp>

double d1(double S, double K, double T, double q, double r, int N, double sigma);
double d2(double S, double K, double T, double q, double r, int N, double sigma);
double normalDist(double d);

double blackScholesCall(double S, double K, double T, double q, double r, int N, double sigma)
{
    // Se^(-qT)N(d1) - Ke^(-rT)*N(d2)
    double d_1= d1(S, K, T, q, r, N, sigma);
    double d_2= d2(S, K, T, q, r, N, sigma);
    return S*std::exp(-q*T)*normalDist(d_1) - K*std::exp(-r*T) * normalDist(d_2);
}

double blackScholesPut(double S, double K, double T, double q, double r, int N, double sigma)
{
    // Se^(-qT)N(d1) - Ke^(-rT)*N(d2)
    double d_1= d1(S, K, T, q, r, N, sigma);
    double d_2= d2(S, K, T, q, r, N, sigma);
    return K*std::exp(-r*T)* normalDist(-d_2) - S*std::exp(-q*T)*normalDist(-d_1);
}


double d1(double S, double K, double T, double q, double r, int N, double sigma)
{
    return (std::log(S/K) + T*(r - q + sigma*sigma/2))/(sigma*std::sqrt(T));
}

double d2(double S, double K, double T, double q, double r, int N, double sigma)
{
    return d1(S, K,T, q, r, N, sigma) - sigma*std::sqrt(T);
}

double normalDist(double d)
{
    boost::math::normal_distribution<double> normalVariable;
    return boost::math::cdf(normalVariable, d);
}
#endif //CPPCODETEST_BLACKSCHOLES_HPP
