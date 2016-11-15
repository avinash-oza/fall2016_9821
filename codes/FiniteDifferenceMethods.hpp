//
// Created by avi on 11/14/16.
//

#ifndef CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
#define CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
#include <Eigen/Dense>

class uFunction
{
    public:
        virtual double evaluate(double t, double x) = 0;

};

class hw8f : public uFunction
{
    public:
        double evaluate(double t, double x)
        {
            return std::exp(x);
        }

};

class hw8gLeft : public uFunction
{
    public:
    double evaluate(double t, double x)
    {
        return std::exp(t - 2.0);
    }

};

class hw8gRight : public uFunction
{
    public:
    double evaluate(double t, double x)
    {
        return std::exp(t + 2.0);
    }

};

//M time intervals, N x intervals

Eigen::VectorXd forwardEuler(uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0, double tFinal, double xLeft, double xRight, int M, int N) {
    // forward euler for the heat equation
    double dx = (xRight - xLeft) / (N);
    double dt = (tFinal - t0) / (M);
    double alpha = (tFinal * N * N) / ((xRight - xLeft) * (xRight - xLeft) * M);
    double F = alpha * dt / (dx * dx);
    double t = t0;
    double x = xLeft;

    VectorXd u(N + 1); // u value at current time
    u.setZero();

    VectorXd u_1(N + 1); // u from the previous time
    u_1.setZero();

    for (int i = 0; i < N + 1; ++i) {
        //set value equal to the x boundary
        u_1[i] = f.evaluate(0, x + i*dx);
        std::cout << "t=" << t << "x=" << x + i*dx << "u=" << u_1[i] << std::endl;
    }

    for (int i = 0; i < M; ++i) {
        t = t0 + i*dt;
        // iterate through each x value
        for (int j = 1; j < N; ++j) {
            u[j] = alpha*u_1[j + 1]  + (1-2*alpha) * u_1[j] + alpha*u_1[j-1];
            std::cout << "t=" << t << " x=" << x + j*dx  << " u= " <<  u[j] << "EXACT:" << std::exp(t + x + j*dx)
                      << "DIFF:" <<  abs(std::exp(t + x + j*dx) - u[j]) << std::endl;
        }



        // set boundary conditions
        u[0] = gLeftFunc.evaluate(t, xLeft);
        u[N] = gRightFunc.evaluate(t, xRight);

        // update u_1 as current values to keep for next iteration
        u_1 = u;
    }
    return u;
}


#endif //CPPCODETEST_FINITEDIFFERENCEMETHODS_HPP
