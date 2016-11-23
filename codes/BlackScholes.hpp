//
// Created by avi on 9/26/16.
//

#ifndef CPPCODETEST_BLACKSCHOLES_HPP
#define CPPCODETEST_BLACKSCHOLES_HPP

#include <cmath>
#include <boost/math/distributions/normal.hpp>

class BlackScholesOption
{
    public:
        BlackScholesOption(double S, double K, double T, double q, double r, double sigma):
        S(S), K(K), T(T), q(q), r(r), sigma(sigma) {};

        double callPrice()
        {
            // Se^(-qT)N(d1) - Ke^(-rT)*N(d2)
            double d_1= d1(S, K, T, q, r,  sigma);
            double d_2= d2(S, K, T, q, r, sigma);
            return S*std::exp(-q*T)*normalDist(d_1) - K*std::exp(-r*T) * normalDist(d_2);
        }

        double putPrice()
        {
            // Se^(-qT)N(d1) - Ke^(-rT)*N(d2)
            double d_1= d1(S, K, T, q, r, sigma);
            double d_2= d2(S, K, T, q, r, sigma);
            return K*std::exp(-r*T)* normalDist(-d_2) - S*std::exp(-q*T)*normalDist(-d_1);
        }


        double d1(double S, double K, double T, double q, double r, double sigma)
        {
            return (std::log(S/K) + T*(r - q + sigma*sigma/2))/(sigma*std::sqrt(T));
        }

        double d2(double S, double K, double T, double q, double r, double sigma)
        {
            return d1(S, K,T, q, r, sigma) - sigma*std::sqrt(T);
        }

        double normalDist(double d)
        {
            boost::math::normal_distribution<double> normalVariable;
            return boost::math::cdf(normalVariable, d);
        }

        double putDelta()
        {
            double d_1 = d1(S, K, T, q, r, sigma);
            return -1*std::exp(-q*T)*normalDist(-d_1);
        }

    double getS() const {
        return S;
    }

    void setS(double S) {
        BlackScholesOption::S = S;
    }

    double getK() const {
        return K;
    }

    void setK(double K) {
        BlackScholesOption::K = K;
    }

    double getT() const {
        return T;
    }

    void setT(double T) {
        BlackScholesOption::T = T;
    }

    double getQ() const {
        return q;
    }

    void setQ(double q) {
        BlackScholesOption::q = q;
    }

    double getR() const {
        return r;
    }

    void setR(double r) {
        BlackScholesOption::r = r;
    }

    double getSigma() const {
        return sigma;
    }

    void setSigma(double sigma) {
        BlackScholesOption::sigma = sigma;
    }


private:
    double S;
    double K;
    double T;
    double q;
    double r;
    double sigma;

};

#endif //CPPCODETEST_BLACKSCHOLES_HPP
