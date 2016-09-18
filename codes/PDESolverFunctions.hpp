//
// Created by avi on 9/18/16.
//

#ifndef CPPCODETEST_PDESOLVERFUNCTIONS_HPP
#define CPPCODETEST_PDESOLVERFUNCTIONS_HPP

#include <tuple>
#include <Eigen/Dense>

using namespace Eigen;

std::tuple<double, double> getStepSizes(double startX, double endX, double startY, double endY, int M, int N) {
    double xStepSize= (endX - startX) / (N + 1);
    double yStepSize= (endY - startY) / (M + 1);// need to add the end points// need to add the end points

    return std::make_tuple(xStepSize, yStepSize);
}

std::tuple<VectorXd, VectorXd> generateMeshPoints(double startX, double endX, double startY, double endY, int M, int N)
{
    double xStepSize;
    double yStepSize;
    std::tie(xStepSize, yStepSize) = getStepSizes(startX, endX, startY, endY, M, N);

    VectorXd xCoordinates(N);
    for (int i=0; i < N; i++)
    {
        xCoordinates(i) = startX +(i + 1)*xStepSize;
    }

    VectorXd yCoordinates(N);
    for (int j=0; j < M; j++)
    {
        yCoordinates(j) = startY +(j + 1)*yStepSize;
    }

    return std::make_tuple(xCoordinates, yCoordinates);
}

double f(double x, double y)
{
    // function f given in hw3
    double temp1 = (x*x + y*y -2)* (std::sin(x)*std::sin(y));
    double temp2 = 2.0*x*std::cos(x)*std::sin(y);
    double temp3 = 2.0*y*std::sin(x)*std::cos(y);

    return temp1 - temp2 - temp3;
}

double u(double x, double y)
{
    // The u(x, y) function given in hw3
    return 0.5*(y*y + x*x)*std::sin(x)*std::sin(y);
}

// Get the b vector
VectorXd get_b(const int N, const VectorXd &xCoordinates, const VectorXd &yCoordinates, double xStepSize, double yStepSize)
{
    // use the xStepSize as the h here
    double h = xStepSize;

    // Generate f_vector to compute bj
    VectorXd f_vec(N);
    MatrixXd b_all(N,N);	// Using matrix to store all N bj

    // Get the 1st column of b_all
    for (int j = 0; j < N; ++j)
    {
        for (int k = 0; k < N; ++k)
        {
            f_vec(k) = f(xCoordinates(k), yCoordinates(0));
        }
        //b_all(j, (N - 1)) = ((N + 1)*(N + 1))*f_vec(j) + u(xCoordinates(j), 1);		// y(N+2)=1
        // correction of b
        b_all(j, 0) = (h*h)*f_vec(j) + u(xCoordinates(j), 0);	// y(1)=0;
        b_all(0, 0) += u(0, yCoordinates(0));
        b_all((N - 1), 0) += u(1, yCoordinates(0));				// x(N+2)=1
    }


    for (int i = 1; i < N-1; ++i)
    {
        // Get the middle columns of b_all
        for (int j = 0; j < N; ++j)
        {
            // Get the f_vector
            for (int k = 0; k < N; ++k)
            {
                f_vec(k) = f(xCoordinates(k), yCoordinates(i));
            }

            //b_all(j, i) = ((N + 1)*(N + 1))*f_vec(j);
            // Correction of b
            b_all(j, i) = (h*h)*f_vec(j);
            b_all(0, i) += u(0, yCoordinates(i));			// x(1)=0
            b_all((N - 1), i) += u(1, yCoordinates(i));		// x(N+2)=1
        }

    }

    // Get the last column of b_all
    for (int j = 0; j < N; ++j)
    {
        for (int k = 0; k < N; ++k)
        {
            f_vec(k) = f(xCoordinates(k), yCoordinates(N - 1));
        }
        //b_all(j, (N - 1)) = ((N + 1)*(N + 1))*f_vec(j) + u(xCoordinates(j), 1);		// y(N+2)=1
        // correction of b
        b_all(j, (N - 1)) = (h*h)*f_vec(j) + u(xCoordinates(j), 1);
        b_all(0, (N - 1)) += u(0, yCoordinates(N - 1));
        b_all((N - 1), (N - 1)) += u(1, yCoordinates(N - 1));
    }


    // Merge the matrix of b_all to a column vector
    VectorXd b(N*N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = N*i; j < N*(i + 1); ++j)
        {
            if (i > 0)
            {
                b(j) = b_all((j % (N*i)), i);
            }
            else
            {
                b(j) = b_all(j, i);		// i=0
            }

        }
    }

    return b;
}


#endif //CPPCODETEST_PDESOLVERFUNCTIONS_HPP
