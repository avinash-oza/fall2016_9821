#include "Eigen/Dense"
#include "MatrixDecomposition.hpp"
#include "IterativeMethods.hpp"
#include "StoppingCriterionSolvers.hpp"
#include "Solvers.hpp"
#include <tuple>
#include <iomanip>
#include <cmath>

using namespace Eigen;
using namespace std;

void verifyCholeskyDecomposition();
MatrixXd generateTestBandedMatrix(int matrixSize); //defined at end

MatrixXd generateHW3Matrix(int N);

MatrixXd Some(const MatrixXd& x)
{
	return x;
};


int main() {
    /// Keep this line to make the decimals always print out
    std::cout << std::fixed << std::setprecision(6);
    //// Example 1
    //// Forward substitution for discounting rates
    //MatrixXd m1(4, 4);
    //m1 <<	100, 0, 0, 0,		6, 106, 0, 0,	 8, 8, 108, 0,		5, 5, 5, 105;
    //VectorXd v1(4);
    //v1 <<	98, 104, 111, 102;
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "|         Forward substitution        |\n";
    //std::cout << "|     Finding discounting factors     |\n";
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "The bond payment matrix is\n"; std::cout << m1 << std::endl << std::endl;
    //std::cout << "The bond prices vector is\n"; std::cout << v1 << std::endl << std::endl;
    //std::cout << "The exponents of the discounting factors are:\n"; std::cout << forward_subst(m1, v1) << std::endl << std::endl;
    //// Remember, the discounting factos are e^(values returned from the function).

    //// Example 2
    //// LU decomposition without pivotting
//	MatrixXd m2(4, 4);
//	m2 << 2,-1,3,0,   -4,5,-7,-2,   -2,10,-4,-7,   4,-14,8,10;
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "|  LU Decomposition Without Pivoting  |\n";
    //std::cout << "|           Full Matrix               |\n";
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "The original matrix is\n"; std::cout << m2 << std::endl << std::endl;
    //std::cout << "The L matrix is:\n" << std::get<0>(lu_no_pivoting(m2)) << std::endl << std::endl;	// L is found in position 0
    //std::cout << "The U matrix is:\n" << std::get<1>(lu_no_pivoting(m2)) << std::endl << std::endl;	// U is found in position 1

    //// Example 3
    //// LU decomposition with row pivoting
    //MatrixXd m3(5, 5);
    //m3 << 1, 2, -7, -1.5, 2, 4, 4, 0, -6, -2, -2, -1, 2, 6, 2.5, 0, -2, 2, -4, 1, 2, 1, 13, -5, 3.5;
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "|   LU Decomposition With Pivoting    |\n";
    //std::cout << "|           Full Matrix               |\n";
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "The original matrix is\n"; std::cout << m3 << std::endl << std::endl;
    //std::cout << "The L matrix is:\n" << std::get<0>(lu_pivoting(m3)) << std::endl << std::endl;	// L is found in position 0
    //std::cout << "The U matrix is:\n" << std::get<1>(lu_pivoting(m3)) << std::endl << std::endl;	// U is found in position 1
    //std::cout << "The U matrix is:\n" << std::get<2>(lu_pivoting(m3)) << std::endl << std::endl;	// P is found in position 2

    //// Example 4
    //// Finding the inverse of a matrix
    //MatrixXd m4(5, 5);
    //m4 << 1, 2, -7, -1.5, 2,	4, 4, 0, -6, -2,	-2, -1, 2, 6, 2.5,		0, -2, 2, -4, 1,	2, 1, 13, -5, 3.5;
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "|         Inverse of a matrix         |\n";
    //std::cout << "|             Full Matrix             |\n";
    //std::cout << "|-------------------------------------|\n";
    //std::cout << "The input matrix is\n" << m4 << endl;
    //std::cout << "The inverse of the matrix is\n" << inverse(m4) << std::endl;
    //std::cout << "The product OF the two matrices is\n" << m4*inverse(m4) << endl;

    // Example 5
    // Finding discounting factors using the LU decomposition
//	MatrixXd m5(4, 4);
//	m5 << 2, 102, 0, 0,		5, 0, 105, 0,		1.5, 1.5, 1.5, 101.5,		0, 6, 0, 106;
//	VectorXd v5(4);
//	v5 << 101.5, 105.5, 101, 106.75;
//	std::cout << "|-------------------------------------|\n";
//	std::cout << "|           LU Decomposition          |\n";
//	std::cout << "|     Finding discounting factors     |\n";
//	std::cout << "|-------------------------------------|\n";
//	std::cout << "The bond payment matrix is\n"; std::cout << m5 << std::endl << std::endl;
//	std::cout << "The bond prices vector is\n"; std::cout << v5 << std::endl << std::endl;
//	std::cout << showpoint << setprecision(8) << "The exponents of the discounting factors are:\n"; std::cout << linear_solve_lu_row_pivoting(m5, v5) << std::endl << std::endl;
    // Remember, the discounting factos are e^(values returned from the function).


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
    GaussSiedelIteration iterationMethod;

//    //coordinates for mesh
    double startX = 0;
    double endX = 1;

    //coordinates for y direction
    double startY = 0;
    double endY = 1;
    int N = 3; // number of points to make on the x axis
    int M = 3; // number of points on the y axis

    double xStepSize = (endX - startX) / N;
    double yStepSize = (endY - startY) / M;

//    VectorXd xCoordinates(N);
//    for (int i=0; i < N; i++)
//    {
//        xCoordinates(i) = startX +i*xStepSize;
//    }
//
//    VectorXd yCoordinates(N-1);
//    for (int j=0; j < M; j++)
//    {
//        yCoordinates(j) = startY +j*yStepSize;
//    }MatrixXd T;
    VectorXd b(9);
    VectorXd x0(9);
    x0.setZero();
    b << 25, 50 ,150 ,0,0, 50, 0,0, 25;
    MatrixXd T = generateHW3Matrix(N);

    std::cout <<"\n" << T << std::endl;
    std::cout << linear_solve_cholesky(T, b) << std::endl;

    cout << "ret_new:\n" << std::endl;
    cout << consecutive_approximation_solver(T, b, x0, std::pow(10, -6), iterationMethod, 1.15);

	return 0;
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


double hw3FiniteDifferenceMethod(double x, double y)
{
    double temp1 = (x*x + y*y -2)* (std::sin(x)*std::sin(y));
    double temp2 = 2.0*x*std::cos(x)*std::sin(y);
    double temp3 = 2.0*y*std::sin(x)*std::cos(y);

    return temp1 - temp2 - temp3;
}

void verifyCholeskyDecomposition()
{
    // Example 6
    // Cholesky decomposition
//    Original matrix:
//    9 -3  6 -3
//    -3  5 -4  7
//    6 -4 21  3
//    -3  7  3 15
//    The Cholesky factor is:
//    3 -1  2 -1
//    0  2 -1  3
//    0  0  4  2
//    0  0  0  1
//    As a check, U^t*U - A = 0:
//    0 0 0 0
//    0 0 0 0
//    0 0 0 0
//    0 0 0 0
    MatrixXd m6(4, 4);
    m6 << 9, -3, 6, -3,
            -3, 5, -4, 7,
            6, -4, 21, 3,
            -3, 7, 3, 15;
    std::cout << "|-------------------------------------|\n";
    std::cout << "|     Cholesky Decomposition          |\n";
    std::cout << "|-------------------------------------|\n";
    std::cout << "Original matrix: " << std::endl;
    std::cout << m6 << std::endl;
    std::cout << "The Cholesky factor is: " << std::endl;
//     Make new variable with cholesky factor
    MatrixXd choleskyFactor(cholesky(m6));
    std::cout << choleskyFactor << std::endl;
    std::cout << "As a check, U^t*U - A = 0:" << std::endl;
    std::cout << choleskyFactor.transpose()*choleskyFactor - m6 << std::endl;

    VectorXd b(4);
    b.setRandom();
    std::cout << "|-------------------------------------|\n";
    std::cout << "|     Cholesky Linear Solver          |\n";
    std::cout << "|-------------------------------------|\n";
    std::cout << "The difference between the original b and the calculated b is:" << std::endl;
    std::cout << m6*linear_solve_cholesky(m6, b) - b << std::endl;

    MatrixXd A = generateTestBandedMatrix(7);

    std::cout << "A is: " << std::endl;
    std::cout << A << std::endl;
    MatrixXd choleskyFactorBanded (cholesky_banded_matrix(A, 2));

    VectorXd b1(7);
    // Make up some b
    b1.setRandom();
    // Check to see that U^t*U - A = 0
    std::cout << (choleskyFactorBanded.transpose())* choleskyFactorBanded - A << std::endl;

    std::cout << "Verify cholesky banded function against cholesky:" << std::endl;
    std::cout << choleskyFactorBanded - cholesky(A) << std::endl;

    std::cout << "Verify banded cholesky solver vs unbanded solver:" << std::endl;
    std::cout << linear_solve_banded_cholesky(A, b1, 2) - linear_solve_cholesky(A, b1);
}

// generates a banded matrix for testing of band size 2
MatrixXd generateTestBandedMatrix(int matrixSize)
{
    MatrixXd A(matrixSize, matrixSize);
    A.setZero();
    //Generate values for the matrix with band size 2
    std::srand(100);
    for (int i = 0; i < matrixSize - 2; i++)
    {
        int num = std::rand() % 25;
        A(i,i + 1) = num;
        A(i + 1, i) = num;

        int num2 = std::rand() % 15;
        A(i,i + 2) = num2;
        A(i + 2, i) = num2;
    }
    // generate the diagnoal entries
    for (int i = 0; i < matrixSize; i++)
    {
        A(i,i) = 20 + std::rand() % 50;
    }

    return A;
}

MatrixXd generateTestIterativeMethodMatrix()
{
    Eigen::MatrixXd A=Eigen::MatrixXd::Zero(14,14);
    Eigen::VectorXd b=Eigen::VectorXd::Zero(14);
    Eigen::VectorXd x_0=Eigen::VectorXd::Zero(14);
    double tol=0.000001;

    for(int j=0; j<14; ++j)
    {
        b(j)=j*j;
        x_0(j)=1;
        for(int k=0; k<14; ++k)
        {
            if(j==k)
                A(j,k)=2;

            else if(j==k+1)
                A(j,k)=-1;
            else if(j==k-1)
                A(j,k)=-1;
            else
                continue;
        }
    }
    return A;
}