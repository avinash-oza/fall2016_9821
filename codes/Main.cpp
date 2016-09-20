#include "Eigen/Dense"
#include "IterativeMethods.hpp"
#include "StoppingCriterionSolvers.hpp"
#include "Solvers.hpp"
#include "Examples.hpp"
#include "PDESolverFunctions.hpp"
#include <tuple>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <vector>

using namespace Eigen;
using namespace std;

#define EIGEN_INITIALIZE_MATRICES_BY_NAN

//defined at end
VectorXd generateBVector(int N, VectorXd xMesh, VectorXd yMesh);
MatrixXd generateHW3Matrix(int N);
double u(double x, double y);
VectorXd calculateuExactVector(int N, const VectorXd &xCoordinates, const VectorXd &yCoordinates);
void runOneIteration(int N, double w);
void printCSVMatrix(std::string stringToPrint, const MatrixXd& myMatrix);

int main() {
    /// Keep this line to make the decimals always print out
    std::cout << std::fixed << std::setprecision(9);
//    decompositionExamples();
//    verifyCholeskyDecomposition();

    // Iterative methods test code (Problem 3 from HW3)
//    Eigen::MatrixXd A = MatrixXd::Zero(14, 14);
//    Eigen::VectorXd b = VectorXd::Zero(14);
//    Eigen::VectorXd x_0 = VectorXd::Zero(14);
//    double tol = 0.000001;
//
//    for (int j = 0; j < 14; ++j) {
//        b(j) = j * j;
//        x_0(j) = 1;
//        for (int k = 0; k < 14; ++k) {
//            if (j == k)
//                A(j, k) = 2;
//
//            else if (j == k + 1)
//                A(j, k) = -1;
//            else if (j == k - 1)
//                A(j, k) = -1;
//            else
//                continue;
//        }
//    }

//    SORIteration jacobiIterator;
//
//    cout << "ret_new:\n" << std::endl;
//    cout << consecutive_approximation_solver(A, b, x_0, std::pow(10, -6), jacobiIterator, 1.15);
//
//    vector<int> NtoRunFor= {2};
//    for (double w = 1.02; w < 1.04; w+=0.02) {
//        cout << w << "|";
//        for (auto it = NtoRunFor.begin(); it != NtoRunFor.end(); it++) {
//            runOneIteration(*it, w);
            cout << ",";
//        }
//        cout << endl;
//    }
//    MatrixXd problem1A = MatrixXd(9,9);
//    for(int i = 0 ; i < 9; i++)
//    {
//        problem1A(i,i) = 2;
//    }
//    for(int i = 1 ; i < 9; i++)
//    {
//        problem1A(i,i-1) = 3;
//    }
//    for(int i = 2 ; i < 9; i++)
//    {
//        problem1A(i,i-2) = 4;
//    }
//    for(int i = 0 ; i < 6; i++)
//    {
//        problem1A(i,i+2) = -1;
//    }
//    std::cout << "Probielmn 1a" << std::endl << problem1A << std::endl;
//
//    MatrixXd L,U,P;
//    std::tie(L,U,P) = lu_pivoting(problem1A);
//    printCSVMatrix("L", L);
//    printCSVMatrix("U", U);
//    printCSVMatrix("P", P);
//
    return 0;
}

void runOneIteration(int N, double w) {
    //Start building the mesh
    double startX = 0;
    double endX = 1;

    //coordinates for y direction
    double startY = 0;
    double endY = 1;

    // h for x and y
    double xStepSize, yStepSize;

    int M = N;

    VectorXd xCoordinates, yCoordinates;

    std::tie(xCoordinates, yCoordinates) = generateMeshPoints(startX, endX, startY, endY, N, M);
    std::tie(xStepSize, yStepSize) = getStepSizes(startX, endX, startY, endY, M, N);

    // Finish building the mesh

    // Construct the interation method and initalize the initial vector
    SORIteration iterationMethod;
    VectorXd x0(N*N);
    x0.setZero();

    // generate the matrix T
    MatrixXd T = generateHW3Matrix(N);
    VectorXd b = get_b(N, xCoordinates, yCoordinates, xStepSize, yStepSize);

    // calculate the exact values we expect
    VectorXd uExact = calculateuExactVector(N, xCoordinates, yCoordinates);

    time_t startTime(time(NULL));

    VectorXd iterativeMethodResult(residual_based_solver(T, b, x0, std::pow(10, -6), iterationMethod, w));

    time_t endTime(time(NULL));

    cout << N << "," << residualError(T, iterativeMethodResult, b);
    cout << maxApproximationError(iterativeMethodResult, uExact);
//    std::cout << "The time elapsed was " << endTime - startTime << " seconds" << std::endl;

}

VectorXd calculateuExactVector(int N, const VectorXd &xCoordinates, const VectorXd &yCoordinates) {
    // calculates the given function u over the x and y coordinates given
    VectorXd uExact(N * N);
    int vectorLocation = 0;

    for (int yIndex = 0 ; yIndex < N; yIndex++)
    {
        for (int xIndex = 0; xIndex < N; xIndex++)
        {
            // MAKE SURE TO CHANGE THIS FUNCTION IF NEEDED
            uExact(vectorLocation) = u(xCoordinates(xIndex), yCoordinates(yIndex));
//            std::cout << "Mesh point is (" << xCoordinates(xIndex) << ", " << yCoordinates(yIndex) << ")" << std::endl;
            vectorLocation++;
        }
    }
    return uExact;
}

MatrixXd generateHW3Matrix(int N) {
    int NSquared = N*N;
    MatrixXd T(NSquared, NSquared);
    T.setZero();

    for (int i = 0; i < N * N; i++) {
        T(i, i) = 4.0;
    }

    for (int j = 1; j < N * N; j++) {
        int valueToCheck = (j + 1)-1;
        if (valueToCheck % N != 0) {
            T(j, j - 1) = -1.0;
        }
    }

    for (int j = 0; j < N * N - 1; j++) {
        int valueToCheck = j + 1;
        if (valueToCheck % N != 0) {
            //assign based on original index, check based on 1 indexed
            T(j, j + 1) = -1.0;
        }
    }

    for (int j = 0; j < N * N - N; j++) {
        T(j, j + N) = -1.0;
    }

    for (int j = N; j < N * N; j++) {
        T(j, j - N) = -1.0;
    }
    return T;

}

void printCSVMatrix(std::string stringToPrint, const MatrixXd& myMatrix)
{
    IOFormat csvFormatter(9, 0, ",");
    std::cout << stringToPrint << std::endl;
    std::cout << myMatrix.format(csvFormatter) << std::endl;
}