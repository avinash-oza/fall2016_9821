//
// Created by avi on 11/14/16.
//

#ifndef CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
#define CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
#include <Eigen/Dense>

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

class hw8gLeft : public uFunction
{
    public:
    double evaluate(double x, double t)
    {
        return std::exp(t - 2.0);
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

//Eigen::VectorXd backwardEuler(uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0, double tFinal, double xLeft, double xRight, int M, int N)
//{
//    double dx = (xRight - xLeft) / (N);
//    double dt = (tFinal - t0) / (M);
//    double alpha = (tFinal * N * N) / ((xRight - xLeft) * (xRight - xLeft) * M);
//    double F = alpha * dt / (dx * dx);
//    double t = t0;
//    double x = xLeft;
//
//    //////
//    VectorXd u(N + 1); // u value at current time
//    u.setZero();
//
//    VectorXd u_1(N + 1); // u from the previous time
//    u_1.setZero();
//    //////
//
//    Eigen::MatrixXd A(N - 1, N - 1);
//    A.setZero();
//
//    for(int i = 0; i < N - 1; ++i)
//    {
//        A(i, i) = 1 + 2*alpha;
//    }
//    for(int i = 0; i < N - 2; ++i)
//    {
//        A(i, i + 1) = - alpha;
//    }
//    for(int i = 1; i < N - 1; ++i)
//    {
//        A(i, i - 1) = - alpha;
//    }
//    std::cout << A << std::endl;
//
//    ///set initial conditions
//    for (int i = 0; i < N + 1; ++i) {
//        //set value equal to the x boundary
//        u_1[i] = f.evaluate(0, x + i*dx);
//        std::cout << "t=" << t << "x=" << x + i*dx << "u=" << u_1[i] << std::endl;
//    }
//    ///
//
///////
//    for (int i = 0; i < M; ++i) {
//        t = t0 + i * dt;
//        // iterate through each x value
//        for (int j = 1; j < N; ++j) {
//            u[j] = alpha * u_1[j + 1] + (1 - 2 * alpha) * u_1[j] + alpha * u_1[j - 1];
//            std::cout << "t=" << t << " x=" << x + j * dx << " u= " << u[j] << "EXACT:" << std::exp(t + x + j * dx)
//                      << "DIFF:" << abs(std::exp(t + x + j * dx) - u[j]) << std::endl;
//        }
//
//
//
//        // set boundary conditions
//        u[0] = gLeftFunc.evaluate(t, xLeft);
//        u[N] = gRightFunc.evaluate(t, xRight);
//
//        // update u_1 as current values to keep for next iteration
//        u_1 = u;
//    }
/////

//    return VectorXd();
//}

//M time intervals, N x intervals

class Mesh
{
public:
    Mesh() {};
    Mesh(double t0, double tFinal, double xLeft, double xRight, int M, int N):
            _t0(t0), _tFinal(tFinal), _xLeft(xLeft), _xRight(xRight), M(M), N(N)
    {
        // construct mesh here
        timeCoordinates = Eigen::VectorXd::Zero(M + 1);
        xCoordinates = Eigen::VectorXd::Zero(N + 1);
        double dx = (_xRight - _xLeft) / N;
        double dt = (_tFinal - _t0) / M;

        for(int i = 0; i < M + 1; ++i)
        {
            timeCoordinates(i) = t0 + i*dt;
        }

        for(int j = 0; j < N + 1; ++j)
        {
            // alter only the x coordinate
            xCoordinates(j) = xLeft + j*dx;
        }

//        std::cout << timeCoordinates << std::endl;
//        std::cout << xCoordinates << std::endl;

    };

    double getX(double index)
    {
        return xCoordinates(index);
    }

    double getT(double index)
    {
        return timeCoordinates(index);
    }
private:
    VectorXd timeCoordinates;
    VectorXd xCoordinates;
    double _t0;
    double _tFinal;
    double _xLeft;
    double _xRight;
    long M;
    long N;

};

class PDESolver
{
public:
    PDESolver(uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0, double tFinal, double xLeft, double xRight, int M, int N) :
            _gLeftFunc(gLeftFunc), _gRightFunc(gRightFunc),
            _f(f), _t0(t0), _tFinal(tFinal), _xLeft(xLeft), _xRight(xRight), M(M), N(N) {
        mesh = Mesh(0, tFinal, xLeft, xRight, M, N);
    };


    Eigen::MatrixXd forwardEuler()
    {
        Eigen::MatrixXd valuesAtNodes = Eigen::MatrixXd::Zero(M + 1, N + 1);

        valuesAtNodes(0,0) = _gLeftFunc.evaluate(_xLeft,0);
        valuesAtNodes(0,N) = _gRightFunc.evaluate(_xRight,0);
        double dx = (_xRight - _xLeft) / N;
        double dt = (_tFinal - _t0) / M;
        double alpha = dt/(dx*dx);

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
private:
    uFunction & _gLeftFunc;
    uFunction & _gRightFunc;
    uFunction & _f;
    Mesh mesh;
    double _t0;
    double _tFinal;
    double _xLeft;
    double _xRight;
    long M;
    long N;


};



#endif //CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
