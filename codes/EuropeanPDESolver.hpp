//
// Created by avi on 11/28/16.
//

#ifndef CPPCODETEST_EUROPEANPDESOLVER_HPP
#define CPPCODETEST_EUROPEANPDESOLVER_HPP

#include <Eigen/Dense>
#include "PDESolver.hpp"
#include "BlackScholes.hpp"

class EuropeanPutPDESolver : public PDESolver {
public:
    EuropeanPutPDESolver(uOptionFunction &gLeftFunc, uOptionFunction &gRightFunc, uOptionFunction &f, double t0,
                         double S0, double K, double T, double q, double r, double sigma, int M, double alphatemp) :
            PDESolver(gLeftFunc, gRightFunc, f, 0, T, 0, 0, M, 0), S0(S0), K(K), T(T), q(q), r(r), sigma(sigma) {
        // set  values based on logic required
        set_xLeft(log(S0 / K) + (r - q - sigma * sigma / 2) * T - 3 * sigma * sqrt(T));
        set_xRight(log(S0 / K) + (r - q - sigma * sigma / 2) * T + 3 * sigma * sqrt(T));
        a = (r-q)/(sigma*sigma) - 0.5;
        b = pow((r-q)/(sigma*sigma) + 0.5, 2) + 2*q/(sigma*sigma);
        _tFinal = T*sigma*sigma/2;

        N = std::floor((_xRight-_xLeft)/sqrt(_tFinal/(alphatemp*M)));
        mesh = Mesh(0, _tFinal, _xLeft, _xRight, M, N);
    };

    double getXCompute()
    {
        return log(S0/K);
    }

    int getxComputeLowerBound() {
        double valueToFind = getXCompute();
        int i = 0;
        for (i; i < N + 1; ++i)
        {
            if (valueToFind < mesh.getX(i))
                break;
        }
        return i - 1;
    }

    double calculateVapprox(MatrixXd &approximations)
    {
        VectorXd V = approximations.row(M);

        double xCompute = getXCompute();
        double lowerxComputePoint = getxComputeLowerBound();

        double xComputeLowerBoundValue = mesh.getX(lowerxComputePoint);
        double xComputeUpperBoundValue = mesh.getX(lowerxComputePoint + 1);

        double S1 = K*exp(xComputeLowerBoundValue);
        double S2 = K * exp(xComputeUpperBoundValue);

        double V1 = V(lowerxComputePoint)*exp(-a*xComputeLowerBoundValue - b*_tFinal);
        double V2 = V(lowerxComputePoint + 1)*exp(-a*xComputeUpperBoundValue - b*_tFinal);

        return ((S2-S0)*V1 + (S0-S1)*V2)/(S2-S1);
    }

    double calculateVapproxLinearInterpolation(MatrixXd & approximations)
    {
        VectorXd V = approximations.row(M);

        double xCompute = getXCompute();
        double lowerxComputePoint = getxComputeLowerBound();

        double xComputeLowerBoundValue = mesh.getX(lowerxComputePoint);
        double xComputeUpperBoundValue = mesh.getX(lowerxComputePoint + 1);

//                ((solver.mesh.getX(i+1)-x)*EulerResultOption(M,i)+(x-solver.mesh.getX(i))*EulerResultOption(M,i+1))/(solver.mesh.getX(i+1)-solver.mesh.getX(i));
        double part1 = (xComputeUpperBoundValue- xCompute)*V(lowerxComputePoint)+(xCompute-xComputeLowerBoundValue)*V(lowerxComputePoint +1);
        double part2 = xComputeUpperBoundValue - xComputeLowerBoundValue;
        return part1/part2;
    }

    VectorXd calculateSpotPrices()
    {
        VectorXd S = VectorXd::Zero(N+1);
        for(int j = 0; j < N+1; ++j)
        {
            S(j) = K*exp(mesh.getX(j));
        }
        return S;
    }

    VectorXd calculateVApproxVector(MatrixXd &approximations)
    {
        VectorXd boundaryApproximations = approximations.row(M); // the lower most row of matrix
        VectorXd S = calculateSpotPrices();
        VectorXd Vappro = VectorXd::Zero(N+1);
        for(int j = 0; j < N+1; ++j)
        {
            S(j) = K*exp(mesh.getX(j));
            Vappro(j) = boundaryApproximations(j)*exp(-a*mesh.getX(j)-b*_tFinal);
        }

        return Vappro;
    }

    double RootMeanSquaredError(MatrixXd& approximations, BlackScholesOption &option)
    {
        double oldSpot = option.getS();
        VectorXd boundaryApproximations = approximations.row(M); // the lower most row of matrix
        long N = getN();
        VectorXd Vappro = calculateVApproxVector(approximations);
        VectorXd S = calculateSpotPrices();
        double Nrms = 0;
        double sum = 0;
        for(int j = 0; j < N+1; ++j)
        {
            option.setS(S(j));
            double price = option.putPrice();
            if(price > 0.00001*S0)
            {
                ++Nrms;
                sum = sum + pow(Vappro(j) - price , 2)/pow(price,2);
            }
        }
        // put the old spot back
        option.setS(oldSpot);
        return sqrt(sum/Nrms);
    }

    double calculateDelta(MatrixXd &approximations)
    {
        int i = getxComputeLowerBound();
        VectorXd Vapprox = calculateVApproxVector(approximations);
        VectorXd S = calculateSpotPrices();

        return (Vapprox(i+1) - Vapprox(i))/(S(i+1) - S(i));
    }

    double calculateGamma(MatrixXd &approximations)
    {
        int i = getxComputeLowerBound();
        VectorXd Vapprox = calculateVApproxVector(approximations);
        VectorXd S = calculateSpotPrices();

        double delta2 = (Vapprox(i+2) - Vapprox(i+1))/(S(i+2) - S(i+1));
        double delta0 = (Vapprox(i) - Vapprox(i-1))/(S(i) - S(i-1));
        return (delta2-delta0)/((S(i+2)+S(i+1)-S(i)-S(i-1))/2);
    }

    double calculateVApprox1(MatrixXd &approximations)
    {
        VectorXd boundaryApproximations = approximations.row(M); // the last row of the matrix
        int i = getxComputeLowerBound();
        double S2 = K*exp(mesh.getX(i+1));
        double S1 = K*exp(mesh.getX(i));
        double V1 = boundaryApproximations(i)*exp(-a*mesh.getX(i)-b*_tFinal);
        double V2 = boundaryApproximations(i+1)*exp(-a*mesh.getX(i+1)-b*_tFinal);

        return ((S2-S0)*V1 + (S0-S1)*V2)/(S2-S1);

    }

    double calculateVApprox2(MatrixXd &approximations)
    {
        double xComputeLowerBoundValue = getXCompute();
        double uApprox = calculateVapproxLinearInterpolation(approximations);
        return std::exp(-a*xComputeLowerBoundValue-b*_tFinal)*uApprox;
    }

    double calculateTheta(MatrixXd &approximations)
    {
        int i = getxComputeLowerBound();
        VectorXd priorboundaryApproximations = approximations.row(M - 1); // the 2nd to last row of the matrix
        VectorXd boundaryApproximations = approximations.row(M); // the 2nd to last row of the matrix
        double S2 = K*exp(mesh.getX(i+1));
        double S1 = K*exp(mesh.getX(i));
        double V1 = boundaryApproximations(i)*exp(-a*mesh.getX(i)-b*_tFinal);
        double V2 = boundaryApproximations(i+1)*exp(-a*mesh.getX(i+1)-b*_tFinal);

        double Vappro1 = calculateVApprox1(approximations);

        double dT = 2*(mesh.getT(M) - mesh.getT(M-1))/(sigma*sigma);
        double V1dT = priorboundaryApproximations(i)*exp(-a*mesh.getX(i)-b*mesh.getT(M-1));
        double V2dT = priorboundaryApproximations(i+1)*exp(-a*mesh.getX(i+1)-b*mesh.getT(M-1));
        double Vt = ((S2-S0)*V1dT + (S0-S1)*V2dT)/(S2-S1);
        return -1 * (Vappro1 - Vt)/dT;
    }

    double calculateErrorPointwise(MatrixXd &approximations, double Vexact)
    {
        return std::abs(calculateVApprox1(approximations) - Vexact);
    }

    double calculateErrorPointwise2(MatrixXd &approximations, double Vexact)
    {
        return std::abs(calculateVApprox2(approximations) - Vexact);
    }


public:
    double S0;
    double q;
    double K;
    double T;
    double r;
    double sigma;
    // constants used in calculation of option price
    double a;
    double b;

};


#endif //CPPCODETEST_EUROPEANPDESOLVER_HPP
