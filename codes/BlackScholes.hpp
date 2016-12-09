//
// Created by avi on 9/26/16.
//

#ifndef CPPCODETEST_BLACKSCHOLES_HPP
#define CPPCODETEST_BLACKSCHOLES_HPP

#include <cmath>
#include <boost/math/distributions/normal.hpp>

boost::math::normal_distribution<> myStandard(0.0, 1.0);	// Standard normal
const long double PI = 3.14159265358979323846;

class BlackScholesOption
{
    public:
        BlackScholesOption(double S, double K, double T, double q, double r, double sigma):
        S(S), K(K), T(T), q(q), r(r), sigma(sigma) {};

        virtual double price() = 0;
        virtual double delta() = 0;

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

    virtual double theta() = 0;


    // Calculating the gamma of the call
    double gamma()
    {

        double temp, d_1;
        d_1 = d1(S, K, T, q, r, sigma);
        return exp(-q*(T))*exp(-0.5*pow(d_1, 2))*1.0 / (S*sigma*sqrt((T) * 2 * PI));;
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


protected:
    double S;
    double K;
    double T;
    double q;
    double r;
    double sigma;

};

class BlackScholesCallOption : public BlackScholesOption
{
public:
    using BlackScholesOption::BlackScholesOption;

    virtual double price() {
        return callPrice();
    }

    virtual double delta() {
        return 0; // needs to be implemented...
    }

    // Calculating the theta of the call
    virtual double theta()
    {
        double temp = sigma*sqrt(T);
        double d_1, d2;
        d_1 = d1(S, K, T, q, r, sigma);
        d2 = d_1 - temp;
        double term1 = r*normalDist(d2)+sigma*pdf(myStandard,d2)/(2.0*sqrt(T));
        double term2 = K*exp(-r*T);

        double term3 = q*S*exp(-q*T)*cdf(myStandard, d_1);
        return term3 - term1 * term2;
    }

private:
    double callPrice()
    {
        // Se^(-qT)N(d1) - Ke^(-rT)*N(d2)
        double d_1= d1(S, K, T, q, r,  sigma);
        double d_2= d2(S, K, T, q, r, sigma);
        return S*std::exp(-q*T)*normalDist(d_1) - K*std::exp(-r*T) * normalDist(d_2);
    }
};

class BlackScholesPutOption : public BlackScholesOption
{
public:
    using BlackScholesOption::BlackScholesOption;

    virtual double price() {
        return putPrice();
    }

    virtual double delta() {
        return putDelta();
    }

    virtual double theta()
    {
        double temp = sigma*sqrt(T);
        double d1, d2;
        d1 = (log(S / K) + (r - q + pow(sigma, 2.0) / 2.0)*T) / temp;
        d2 = d1 - temp;
        double term1 = -exp(-q*T)*S*sigma*pdf(myStandard, d1) / (2.0*sqrt(T));
        double term2 = r*K*exp(-r*T)*cdf(myStandard, -d2);
        double term3 = q*S*exp(-q*T)*cdf(myStandard, -d1);
        return term1 + term2 + term3;
    }




private:

    double putPrice()
    {
        // Se^(-qT)N(d1) - Ke^(-rT)*N(d2)
        double d_1= d1(S, K, T, q, r, sigma);
        double d_2= d2(S, K, T, q, r, sigma);
        return K*std::exp(-r*T)* normalDist(-d_2) - S*std::exp(-q*T)*normalDist(-d_1);
    }


    double putDelta()
    {
        double d_1 = d1(S, K, T, q, r, sigma);
        return -1*std::exp(-q*T)*normalDist(-d_1);
    }
};

#endif //CPPCODETEST_BLACKSCHOLES_HPP
