//
// Created by avi on 12/3/16.
//

#ifndef CPPCODETEST_BARRIEROPTION_HPP
#define CPPCODETEST_BARRIEROPTION_HPP

#include "BlackScholes.hpp"

class BarrierOption
{
private:
    double S0, K, B, r, sigma, q,T;

public:
    BarrierOption(double Spot, double Strike, double Maturity, double Dividend, double Interest, double Volatility, double Barrier)
            :S0(Spot), K(Strike), B(Barrier), r(Interest), sigma(Volatility), q(Dividend), T(Maturity) {;}

    double Price() const
    {
        double alpha(0.0);
        double result;
        // make sure to change option type here
        BlackScholesCallOption Option1(S0, K, T, q, r, sigma);
        BlackScholesCallOption Option2(B*B / S0, K, T, q, r, sigma);

        alpha = (r - q) / (sigma*sigma) - 0.5;
        result = Option1.price() - pow(B / S0, 2 * alpha)*Option2.price();
        return result;
    }
};


#endif //CPPCODETEST_BARRIEROPTION_HPP
