//
// Created by avi on 11/23/16.
//

#ifndef CPPCODETEST_UFUNCTIONS_HPP
#define CPPCODETEST_UFUNCTIONS_HPP

class uFunction
{
    public:
        virtual double evaluate(double x, double t) = 0;

};

class hw8f : public uFunction
{
    public:
        double evaluate(double x, double t)
        {
            return std::exp(x);
        }

};

class hw8fOption : public uFunction
{
	public:
		double sigma;
		double S0;
		double q;
		double K;
		double T;
		double r;
		double a;
		hw8fOption(double _sigma, double _S0, double _q, double _K, double _T, double _r):
			sigma(_sigma), S0(_S0), q(_q), K(_K), T(_T), r(_r)
		{
			a = (r-q)/(sigma*sigma) - 0.5;
		}
		double evaluate(double x, double t)
		{
			return K * exp(a * x) * std::max(1 - exp(x), 0.0);
		}

};

class hw8gLeft : public uFunction
{
    public:
    double evaluate(double x, double t)
    {
        return std::exp(t - 2.0);
    }

};

class hw8gLeftOption : public uFunction
{
	public:
		double sigma;
		double S0;
		double q;
		double K;
		double T;
		double r;
		double a;
		double b;
		hw8gLeftOption(double _sigma, double _S0, double _q, double _K, double _T, double _r):
			sigma(_sigma), S0(_S0), q(_q), K(_K), T(_T), r(_r)
		{
			a = (r-q)/(sigma*sigma) - 0.5;
			b = pow((r-q)/(sigma*sigma) + 0.5, 2) + 2*q/(sigma*sigma);
		}
		double evaluate(double x, double t)
		{
			return K * exp(a * x + b * t) * (exp(-(2 * r * t) / (sigma * sigma)) - exp(x - 2 * q * t / (sigma * sigma)));
		}
};

class uExact : public uFunction
{
public:
    double evaluate(double x, double t)
    {
        return std::exp(t + x);
    }

};

class hw8gRight : public uFunction
{
    public:
    double evaluate(double x, double t)
    {
        return std::exp(t + 2.0);
    }

};

class hw8gRightOption : public uFunction
{
	public:
    double evaluate(double x, double t)
    {
        return 0.0;
    }
};

#include <Eigen/Dense>

#endif //CPPCODETEST_UFUNCTIONS_HPP
