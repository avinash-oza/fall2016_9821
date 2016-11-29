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



    virtual Eigen::MatrixXd forwardEuler()
    {
        Eigen::MatrixXd valuesAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);

        valuesAtNodes(0,0) = _gLeftFunc.evaluate(_xLeft,0);
        valuesAtNodes(0,N) = _gRightFunc.evaluate(_xRight,0);
        double alpha = getAlpha();

        // create initial vector
        long size = N - 1;
        double deltaX = (_xRight - _xLeft) / N;

        VectorXd U = generateInitalUMatrix(size);
        MatrixXd A = generateForwardEulerMatrix(size);

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
                  valuesAtNodes(timeIndex, i) = calculateUOnMesh(timeIndex, i, U(i-1));
            }
        }
        return valuesAtNodes;
    }



    MatrixXd backwardEuler(LinearSolveMethod linearSolverMethod, double tol, double omega)
    {
        MatrixXd valuesAtNodes = MatrixXd::Zero(M + 1, N + 1);

        long numberTimeSteps = M;

        //create U initial
        long size = N - 1;
        double deltaX = (_xRight - _xLeft) / N;

        VectorXd U = generateInitalUMatrix(size);

        //create A
        double alpha = getAlpha();

        MatrixXd A = generateBackwardEulerMatrix(size);
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

        VectorXd U = generateInitalUMatrix(size);

        for (long i = 1; i < N; ++i)
        {
            valuesAtNodes(0,i) = U(i-1);
        }

        MatrixXd A = generateCrankNicolsonMatrix(size);
        MatrixXd B = getCrankNicolsonBMatrix(size);

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

    virtual double calculateUOnMesh(long timeIndex, long currentIndex,  double europeanUValue)
    {
        // calculates the value for U at a given point on the mesh. Used to handle the calculation of the early exercise premium
        return europeanUValue;
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

    double getAlpha() const {
        double dx = (_xRight - _xLeft) / N;
        double dt = (_tFinal - _t0) / M;
        double alpha = dt/(dx*dx);
        return alpha;
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

    double getM() const {return M;};
    double getN() const {return N;};


protected:
    // to organize the methods better

    MatrixXd generateForwardEulerMatrix(long size) const {
        double alpha = getAlpha();
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
        return A;
    }

    VectorXd generateInitalUMatrix(long size) const {
        VectorXd U = VectorXd::Zero(size);

        for (long i = 0; i < size; ++i)
        {
            U(i) = _f.evaluate(mesh.getX(i + 1), mesh.getT(0));
        }
        return U;
    }

    MatrixXd generateBackwardEulerMatrix(long size) const {
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
        return A;
    }

    MatrixXd generateCrankNicolsonMatrix(long size) const {

        double alpha = getAlpha();
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
        return A;
    }

    MatrixXd getCrankNicolsonBMatrix(long size) const {

        double alpha = getAlpha();
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
        return B;
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


#endif //CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
