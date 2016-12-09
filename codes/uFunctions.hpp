//
// Created by avi on 11/23/16.
//

#ifndef CPPCODETEST_UFUNCTIONS_HPP
#define CPPCODETEST_UFUNCTIONS_HPP

#include <Eigen/Dense>

/**
 * Base function for defining boundary conditions. All boundary functions should inherit from here somehow.
 */
class uFunction
{
    public:
    /**
     * Main entry point for uFunction. This will be the method called for calculation
     * @param x
     * @param t
     * @return
     */
        virtual double evaluate(double x, double t) = 0;
    virtual void setUp()
    {
        if (!isSetUp)
        {
            _setUp();
            isSetUp = true;
        }

    }

    virtual void _setUp() {};

protected:
    bool isSetUp = false;

};
/**
 * uOptionFunction provides an interface for defining boundary conditions driven by option inputs.
 */
class uOptionFunction : public uFunction

{
public:
    uOptionFunction(double _sigma, double _S0, double _q, double _K, double _T, double _r):
            sigma(_sigma), S0(_S0), q(_q), K(_K), T(_T), r(_r)
    { }

    virtual void _setUp()
    {
        seta();
        setb();

    }

    double evaluate(double x, double t)
    {
        return K * exp(a * x) * std::max(1 - exp(x), 0.0);
    }

    virtual void seta()
    {
        a = (r-q)/(sigma*sigma) - 0.5;
    }

    virtual void setb()
    {
        b = pow((r-q)/(sigma*sigma) + 0.5, 2) + 2*q/(sigma*sigma);
    }

    double sigma;
    double S0;
    double q;
    double K;
    double T;
    double r;
    double a;
    double b;
};


/////////////////////////////////////// IMPLEMENTATIONS//////////////////////////////////////////

class hw8f : public uFunction
{
    public:

        double evaluate(double x, double t)
        {
            return std::exp(x);
        }

};

class fOption : public uOptionFunction
{
	public:
        using uOptionFunction::uOptionFunction;
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

class gEuropeanLeft : public uOptionFunction
{
	public:
        using uOptionFunction::uOptionFunction;
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

class gAmericanRight : public uOptionFunction
{
	public:
        using uOptionFunction::uOptionFunction;
        double evaluate(double x, double t)
        {
            return 0.0;
        }
};

class gEuropeanRight : public uOptionFunction
{
public:
    using uOptionFunction::uOptionFunction;
};

class gAmericanLeftFunc : public uOptionFunction
{
public:
    using uOptionFunction::uOptionFunction;
    double evaluate(double x, double t)
    {
        return K*std::exp(a*x + b*t) * (1.0 - std::exp(x));
//        return K*std::exp(a*x + b*t) * (std::exp(x) - 1.0); // CALL
    }
};


#endif //CPPCODETEST_UFUNCTIONS_HPP
