//
// Created by avi on 11/14/16.
//

#ifndef CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
#define CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
#include <Eigen/Dense>
#include "PDEMesh.hpp"
#include "uFunctions.hpp"


using namespace Eigen;
enum LinearSolveMethod {LU, SOR}; // define enum for linear solve

//M time intervals, N x intervals

class PDESolver
{
public:
    PDESolver(uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0, double tFinal, double xLeft, double xRight, int M, int N) :
            _gLeftFunc(gLeftFunc), _gRightFunc(gRightFunc),
            _f(f), _t0(t0), _tFinal(tFinal), _xLeft(xLeft), _xRight(xRight), M(M), N(N) {
        mesh = Mesh(0, tFinal, xLeft, xRight, M, N);
    };



	double getM() const {return M;};
	double getN() const {return N;};
    Eigen::MatrixXd forwardEuler()
    {
        Eigen::MatrixXd valuesAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);

        valuesAtNodes(0,0) = _gLeftFunc.evaluate(_xLeft,0);
        valuesAtNodes(0,N) = _gRightFunc.evaluate(_xRight,0);
        double alpha = getAlpha();

        // create initial vector
        long size = N - 1;
        double deltaX = (_xRight - _xLeft) / N;

        VectorXd U = VectorXd::Zero(size);

        for (long i = 0; i < size; ++i)
        {
            U(i) = _f.evaluate(mesh.getX(i + 1), mesh.getT(0));
        }

        //create constant Euler matrix A

        MatrixXd A = MatrixXd::Zero(size, size);

        for (long i = 0; i < size; ++i)
        {
            A(i,i) = 1 - 2 * alpha;
        }

        for (long i = 0; i < size - 1; ++i)
        {
            A(i, i + 1) = alpha;
        }

        for (long i = 1; i < size; ++i)
        {
            A(i, i - 1) = alpha;
        }
        //

        for (long i = 1; i < N; ++i)
        {
            valuesAtNodes(0,i) = U(i-1);
        }

        for (long timeIndex = 1; timeIndex <= M; ++timeIndex)
        {
            // create b vector

            VectorXd b = VectorXd::Zero(size);

            b(0) = alpha * _gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex - 1));
            b(size - 1) = alpha * _gRightFunc.evaluate(_xRight, mesh.getT(timeIndex - 1));

            U = A * U + b;
            //
            valuesAtNodes(timeIndex,0) = _gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex));
            valuesAtNodes(timeIndex,N) = _gRightFunc.evaluate(_xRight, mesh.getT(timeIndex));

            for (long i = 1; i < N; ++i)
            {
                valuesAtNodes(timeIndex,i) = U(i-1);
            }
        }
        return valuesAtNodes;
    }

    double getAlpha() const {
        double dx = (_xRight - _xLeft) / N;
        double dt = (_tFinal - _t0) / M;
        double alpha = dt/(dx*dx);
        return alpha;
    }

    VectorXd CrankNicolsonCreateBvector(long timeIndex)
    {
        long size = N - 1;
        double alpha = getAlpha();

        VectorXd b = VectorXd::Zero(size);

        b(0) = alpha / 2.0 * (_gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex)) + _gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex - 1)));
        b(size - 1) = alpha / 2.0 * (_gRightFunc.evaluate(_xRight, mesh.getT(timeIndex)) + _gRightFunc.evaluate(_xRight,mesh.getT(timeIndex - 1)));

        return b;
    }


    double RootMeanSquaredError(MatrixXd& approximations, uFunction &uExact)
    {
        double dx = (_xRight - _xLeft) / N;
        double uExactValue;

        VectorXd boundaryApproximations = approximations.row(M); // the lower most row of matrix
        VectorXd totalScaledError(N + 1); // contains the difference of the error squared divided by the exact value squared

        for (long i = 0; i < N + 1; ++i) {
            uExactValue = uExact.evaluate(mesh.getX(i), _tFinal);
            totalScaledError(i) = std::pow(std::abs((boundaryApproximations(i) - uExactValue)),2) / (std::pow(std::abs(uExactValue), 2));
        }

        return std::sqrt(totalScaledError.sum()/(N + 1));
    }

    MatrixXd backwardEuler(LinearSolveMethod linearSolverMethod, double tol, double omega)
    {

        MatrixXd valuesAtNodes = MatrixXd::Zero(M + 1, N + 1);

        long numberTimeSteps = M;

        //create U initial
        long size = N - 1;
        double deltaX = (_xRight - _xLeft) / N;

        VectorXd U = VectorXd::Zero(size);

        for (long i = 0; i < size; ++i)
        {
            U(i) = _f.evaluate(mesh.getX(i + 1), mesh.getT(0));
        }

        //

        //create A
        double alpha = getAlpha();

        MatrixXd A = MatrixXd::Zero(size, size);

        for (long i = 0; i < size; ++i)
        {
            A(i,i) = 1 + 2 * alpha;
        }

        for (long i = 0; i < size - 1; ++i)
        {
            A(i, i + 1) = -alpha;
        }

        for (long i = 1; i < size; ++i)
        {
            A(i, i - 1) = -alpha;
        }

        //
        VectorXd x0 = VectorXd::Zero(N - 1);

        valuesAtNodes(0,0) = _gLeftFunc.evaluate(_xLeft,0);
        valuesAtNodes(0,N) = _gRightFunc.evaluate(_xRight,0);

        for (long i = 1; i < N; ++i)
        {
            valuesAtNodes(0,i) = U(i-1);
        }

        for (long timeIndex = 1; timeIndex <= numberTimeSteps; ++timeIndex)
        {
            VectorXd b = VectorXd::Zero(size);

            b(0) = alpha * _gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex));
            b(size - 1) = alpha * _gRightFunc.evaluate(_xRight, mesh.getT(timeIndex));


            //
            int ic;

            switch (linearSolverMethod)
            {
                case LinearSolveMethod ::LU :
                    U = linear_solve_lu_no_pivoting(A, U + b);
                    break;
                case LinearSolveMethod ::SOR :
                    SORIteration iterationMethod;
                    std::tie(U, ic) = consecutive_approximation_solver(A, U + b, U, tol, iterationMethod, omega);
                    break;
            }

            valuesAtNodes(timeIndex,0) = _gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex));
            valuesAtNodes(timeIndex,N) = _gRightFunc.evaluate(_xRight, mesh.getT(timeIndex));

            for (long i = 1; i < N; ++i)
            {
                valuesAtNodes(timeIndex,i) = U(i-1);
            }
        }
        return valuesAtNodes;
    }


    Eigen::MatrixXd CrankNicolson(LinearSolveMethod linearSolverMethod, double tol, double omega)
    {
        Eigen::MatrixXd valuesAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);

        valuesAtNodes(0,0) = _gLeftFunc.evaluate(_xLeft,0);
        valuesAtNodes(0,N) = _gRightFunc.evaluate(_xRight,0);
        double alpha = getAlpha();

        // create initial vector
        long size = N - 1;
        double deltaX = (_xRight - _xLeft) / N;
        long numberTimeSteps = M;
        VectorXd x0 = VectorXd::Zero(size);

        VectorXd U = VectorXd::Zero(size);

        for (long i = 0; i < size; ++i)
        {
            U(i) = _f.evaluate(mesh.getX(i + 1), mesh.getT(0));
        }

        //

        for (long i = 1; i < N; ++i)
        {
            valuesAtNodes(0,i) = U(i-1);
        }

        // A

        MatrixXd A = MatrixXd::Zero(size, size);

        for (long i = 0; i < size; ++i)
        {
            A(i,i) = 1 + alpha;
        }

        for (long i = 0; i < size - 1; ++i)
        {
            A(i, i + 1) = -alpha / 2.0;
        }

        for (long i = 1; i < size; ++i)
        {
            A(i, i - 1) = -alpha / 2.0;
        }

        // end A


        // create B

        MatrixXd B = MatrixXd::Zero(size, size);

        for (long i = 0; i < size; ++i)
        {
            B(i,i) = 1 - alpha;
        }

        for (long i = 0; i < size - 1; ++i)
        {
            B(i, i + 1) = alpha / 2;
        }

        for (long i = 1; i < size; ++i)
        {
            B(i, i - 1) = alpha / 2;
        }
        // end B

        for (long timeIndex = 1; timeIndex <= numberTimeSteps; ++timeIndex)
        {
            MatrixXd b = CrankNicolsonCreateBvector(timeIndex);
            int ic;

            switch (linearSolverMethod)
            {
                case LinearSolveMethod ::LU :
                    U = linear_solve_lu_no_pivoting(A, B * U + b);
                    break;
                case LinearSolveMethod ::SOR :
                    SORIteration iterationMethod;
                    std::tie(U, ic) = consecutive_approximation_solver(A, B * U + b, U, tol, iterationMethod, omega);
                    break;
            }

            valuesAtNodes(timeIndex,0) = _gLeftFunc.evaluate(_xLeft, mesh.getT(timeIndex));
            valuesAtNodes(timeIndex,N) = _gRightFunc.evaluate(_xRight,mesh.getT(timeIndex));

            for (long i = 1; i < N; ++i)
            {
                valuesAtNodes(timeIndex,i) = U(i-1);
            }
        }

        return valuesAtNodes;

    }


    double MaxPointwiseApproximationError(MatrixXd &approximations, uFunction& uExact) const {
        double dx = (_xRight - _xLeft)/ N;

        VectorXd boundaryApproximations = approximations.row(M); // the lower most row of matrix
        VectorXd uExactBoundary(N + 1);

        for (long i = 0; i < N + 1; ++i)
        {
            double x = _xLeft + i * dx;
            uExactBoundary(i)  = uExact.evaluate(x, _tFinal);
        }

        VectorXd difference = (boundaryApproximations - uExactBoundary);

        return difference.cwiseAbs().maxCoeff();
    }

    double get_xLeft() const {
        return _xLeft;
    }

    void set_xLeft(double _xLeft) {
        PDESolver::_xLeft = _xLeft;
    }

    double get_xRight() const {
        return _xRight;
    }

    void set_xRight(double _xRight) {
        PDESolver::_xRight = _xRight;
    }

public:
    uFunction & _gLeftFunc;
    uFunction & _gRightFunc;
    uFunction & _f;
    double _t0;
    double _tFinal;
    double _xLeft;
    double _xRight;
    long M;
    long N;
	Mesh mesh;


};

class EuropeanPutPDESolver : public PDESolver {
public:
    EuropeanPutPDESolver(uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0,
                         double S0, double K, double T, double q, double r, double sigma, int M, double alphatemp, BlackScholesOption blackScholesOption1) :
            PDESolver(gLeftFunc, gRightFunc, f, 0, T, 0, 0, M, 0), S0(S0), K(K), T(T), q(q), r(r), sigma(sigma), blackScholesOption(blackScholesOption1) {
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

    double calculateVapprox(double S0, MatrixXd & approximations)
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

    double RootMeanSquaredError(MatrixXd& approximations)
    {
        double oldSpot = blackScholesOption.getS();
        VectorXd boundaryApproximations = approximations.row(M); // the lower most row of matrix
        long N = getN();
        VectorXd Vappro = calculateVApproxVector(approximations);
        VectorXd S = calculateSpotPrices();
        double Nrms = 0;
        double sum = 0;
        for(int j = 0; j < N+1; ++j)
        {
            blackScholesOption.setS(S(j));
            double price = blackScholesOption.putPrice();
            if(price > 0.00001*S0)
            {
                ++Nrms;
                sum = sum + pow(Vappro(j) - price , 2)/pow(price,2);
            }
        }
        // put the old spot back
        blackScholesOption.setS(oldSpot);
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

    double calculateVApprox1(VectorXd &boundaryApproximations)
    {
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

        double Vappro1 = calculateVApprox1(boundaryApproximations);

        double dT = 2*(mesh.getT(M) - mesh.getT(M-1))/(sigma*sigma);
		double V1dT = priorboundaryApproximations(i)*exp(-a*mesh.getX(i)-b*mesh.getT(M-1));
		double V2dT = priorboundaryApproximations(i+1)*exp(-a*mesh.getX(i+1)-b*mesh.getT(M-1));
		double Vt = ((S2-S0)*V1dT + (S0-S1)*V2dT)/(S2-S1);
		return (Vappro1 - Vt)/dT;
    }

    double calculateErrorPointwise(MatrixXd &approximations)
    {
        blackScholesOption.setS(S0);
        VectorXd boundaryApproximations = approximations.row(M);
        return std::abs(calculateVApprox1(boundaryApproximations) - blackScholesOption.putPrice());
    }

    double calculateErrorPointwise2(MatrixXd &approximations)
    {
        blackScholesOption.setS(S0);
        return std::abs(calculateVApprox2(approximations) - blackScholesOption.putPrice());
    }



public:
    double S0;
    double q;
    double K;
    double T;
    double r;
    double sigma;
    BlackScholesOption blackScholesOption;
private:
    // constants used in calculation of option price
    double a;
    double b;

};

#endif //CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
