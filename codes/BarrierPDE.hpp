#ifndef BarrierPDE_HPP
#define BarrierPDE_HPP

#include "BlackScholes.hpp"
#include "PDESolver.hpp"
#include "EuropeanPDESolver.hpp"

#include <cmath>

using namespace std;

class DownOutCallPDESolver :public EuropeanPDESolver
{
public:
    DownOutCallPDESolver(uBarrierOption &LeftFunction, uBarrierOption &RightFunction, uBarrierOption &f, double t0,
                             double Spot, double Strike, double Maturity, double Dividend, double Interest,
                             double Volatility, double AlphaTemp, double M, double Barrier)
            :EuropeanPDESolver(LeftFunction, RightFunction, f, t0, Spot, Strike, Maturity, Dividend, Interest, Volatility, M, AlphaTemp), B(Barrier) {
        ;
    }
    virtual void _setUp()
    {
        set_xLeft(calculateXLeft());
        a = (r - q) / (sigma*sigma) - 0.5;
        b = pow((r - q) / (sigma*sigma) + 0.5, 2.0) + 2 * q / (sigma*sigma);
        _tFinal = T*sigma*sigma / 2.0;

        delta_tau = _tFinal / (M + 0.0);
        double delta_xTemp = sqrt(delta_tau / alphatemp);

        N_left = floor((getXCompute()-get_xLeft())/ delta_xTemp);

        delta_x = (getXCompute() - get_xLeft()) / N_left;

        double x_right_temp = log(S0 / K) + (r - q - sigma*sigma / 2.0)*T + 3 * sigma*sqrt(T);
        double xr = log(S0 / K) + (r - q - sigma*sigma / 2.0)*T + 3 * sigma*sqrt(T);
        int N_right = ceil((x_right_temp - getXCompute()) / delta_x);

        N = N_left + N_right;
        set_xRight(calculateXRight(N_right));

        mesh = Mesh(0, _tFinal, _xLeft, _xRight, M, N);
    }

    virtual double calculateXRight(int N_right) { return getXCompute() + N_right * delta_x; }

    virtual int getNRight() {
        double x_right_temp = log(S0 / K) + (r - q - sigma * sigma / 2.0) * T + 3 * sigma * sqrt(T);
        double xr = log(S0 / K) + (r - q - sigma * sigma / 2.0) * T + 3 * sigma * sqrt(T);
        int N_right = ceil((x_right_temp - getXCompute()) / delta_x);
        return N_right;
    }

    virtual double calculateXLeft() const { return log(B / K); }

    // Approximate function
    double calculateVApprox1(MatrixXd &approximations)
    {
        VectorXd boundaryApproximations = approximations.row(M); // the last row of the matrix
        int xComputeLocation = getxComputeLowerBound();
        double xCompute = getXCompute();
        double V1 = boundaryApproximations(xComputeLocation)*exp(-a*xCompute - b*_tFinal);
        return V1;
    }

    // Delta of the function
    double calculateDelta(MatrixXd &approximations)
    {
        VectorXd boundaryApproximations = approximations.row(M); // the last row of the matrix
        double xCompute = getXCompute();

        // The x's at the provided indices
        double x_minus = mesh.getX(N_left - 1);
        double x_zero = mesh.getX(N_left);
        double x_plus = mesh.getX(N_left + 1);

        // The stock prices with those x's
        double S_minus = K*exp(x_minus);
        double S_zero = K*exp(x_zero);
        double S_plus = K*exp(x_plus);

        // The payoff functions
        double V_minus = exp(-a*(x_minus)-b*_tFinal)* boundaryApproximations(N_left - 1);
        double V_zero = exp(-a*(x_zero)-b*_tFinal)* boundaryApproximations(N_left);
        double V_plus = exp(-a*(x_plus)-b*_tFinal)* boundaryApproximations(N_left + 1);

        return (V_plus - V_minus) / (S_plus - S_minus);
    }

    // Gamma of the function
    double calculateGamma(MatrixXd &approximations)
    {
        VectorXd boundaryApproximations = approximations.row(M); // the last row of the matrix
        double xCompute = getXCompute();

        // The x's at the provided indices
        double x_minus = mesh.getX(N_left - 1);
        double x_zero = mesh.getX(N_left);
        double x_plus = mesh.getX(N_left + 1);

        // The stock prices with those x's
        double S_minus = K*exp(x_minus);
        double S_zero = K*exp(x_zero);
        double S_plus = K*exp(x_plus);

        // The payoff functions
        double V_minus = exp(-a*(x_minus)-b*_tFinal)* boundaryApproximations(N_left - 1);
        double V_zero = exp(-a*(x_zero)-b*_tFinal)* boundaryApproximations(N_left);
        double V_plus = exp(-a*(x_plus)-b*_tFinal)* boundaryApproximations(N_left + 1);

        double N = (S_zero - S_minus)*V_plus - (S_plus - S_minus)*V_zero + (S_plus - S_zero)*V_minus;
        double D = (S_zero - S_minus)*(S_plus - S_zero)*(0.5*(S_plus - S_minus));
        return N / D;
    }

    // Theta of the function
    double calculateTheta(MatrixXd &approximations)
    {
        VectorXd boundaryApproximations = approximations.row(M); // the last row of the matrix
        VectorXd priorBoundaryApproximations = approximations.row(M-1); // the second to last row of the matrix
        double V_zero, V_delta;
        // The x's at the provided indices
        double x_zero = mesh.getX(N_left);

        V_delta = exp(-a*x_zero - b*(_tFinal - delta_tau))*priorBoundaryApproximations(N_left);
        V_zero = calculateVApprox1(approximations);

        return -(V_zero - V_delta) / get_dt();
    }


    int getxComputeLowerBound() {
        return N_left;
    }


    double get_Alpha() const
    {
        return delta_tau / (delta_x*delta_x);
    }

    double get_N() const
    {
        return N;
    }

    double get_delta_x()const
    {
        return delta_x;
    }

    double get_delta_tau() const
    {
        return delta_tau;
    }

    long get_N_left() const
    {
        return N_left;
    }

    double get_dt() const {
        return 2.0*delta_tau/(sigma*sigma);
    }
public:
    double B;

private:
    double a;
    double b;

    // Constants to be reported for homework 9
    double delta_x;
    double delta_tau;
    double alpha;
    long N_left;
};

class DoubleBarrierPDESolver : public DownOutCallPDESolver
{
public:
    DoubleBarrierPDESolver(uBarrierOption &LeftFunction, uBarrierOption &RightFunction, uBarrierOption &f, double t0,
                           double Spot, double Strike, double Maturity, double Dividend, double Interest,
                           double Volatility, double AlphaTemp, double M, double Barrier, double upperBarrier) : DownOutCallPDESolver(
            LeftFunction, RightFunction, f, t0, Spot, Strike, Maturity, Dividend, Interest, Volatility, AlphaTemp, M,
            Barrier), upperBarrier(upperBarrier) {}

    double calculateXRight(int N_right) override {
        return log(upperBarrier/K);
    }

protected:
    double upperBarrier;

};

#endif