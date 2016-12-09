//
// Created by avi on 11/16/16.
//
#include <Eigen/Dense>

#ifndef CPPCODETEST_PDEMESH_HPP
#define CPPCODETEST_PDEMESH_HPP
using namespace Eigen;

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
    };

    double getX(double index) const
    {
        return xCoordinates(index);
    }

    double getT(double index) const
    {
        return timeCoordinates(index);
    }

    VectorXd getTimeCoordinates() const
    {
        return timeCoordinates;
    }

    VectorXd getxCoordinates() const
    {
        return xCoordinates;
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


#endif //CPPCODETEST_PDEMESH_HPP
