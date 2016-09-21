#ifndef Solvers_HPP
#define Solvers_HPP

#include "MatrixDecomposition.hpp"
#include "Eigen/Dense"
#include <tuple>
#include <iostream>

using namespace Eigen;

void printCSVMatrix(std::string stringToPrint, const MatrixXd& myMatrix)
{
	IOFormat csvFormatter(9, 0, ",");
	std::cout << stringToPrint << std::endl;
	std::cout << myMatrix.format(csvFormatter) << std::endl;
}

// Function 1
// Forward substitution 
// Operation count = n^2 + O(n)
VectorXd forward_subst(const MatrixXd& L, const VectorXd b)
{
	int size = L.cols();
	VectorXd x(size);
	x(0) = b(0) / L(0, 0);
	for (int j = 1; j <= size - 1; ++j)
	{
		double sum = 0;
		for (int k = 0; k <= j - 1; ++k)
		{
			sum += L(j, k)*x(k);
		}
		x(j) = (b(j) - sum) / L(j, j);
	}
	return x;
}

// Function 2
// Forward substitution for bidiagonal matrices
// Operation count = 3*n - 2 
VectorXd forward_subst_bidiag(const MatrixXd& L, const VectorXd b)
{
	int size = L.cols();
	VectorXd x(size);
	x(0) = b(0) / L(0, 0);
	for (int j = 1; j <= size - 1; ++j)
	{
		x(j) = (b(j) - L(j,j-1)*x(j-1)) / L(j, j);
	}
	return x;
}

// Function 3
// Backward substitution
// Operation count = n^2 + O(n)
VectorXd backward_subst(const MatrixXd& U, const VectorXd b)
{
	int size = U.cols();
	VectorXd x(size);
	x(size-1) = b(size - 1) / U(size - 1, size - 1);
	for (int j = size-2; j >=0; --j)
	{
		double sum = 0;
		for (int k = j+1; k <= size - 1; ++k)
		{
			sum += U(j, k)*x(k);
		}
		x(j) = (b(j) - sum) / U(j, j);
	}
	return x;
}

// Function 4
// Backward substitution for bidiagonal matrices
// Operation count = 3*n - 2 
VectorXd backward_subst_bidiag(const MatrixXd& U, const VectorXd b)
{
	int size = U.cols();
	VectorXd x(size);
	x(size - 1) = b(size - 1) / U(size - 1, size - 1);
	for (int j = size - 2; j >= 0; --j)
	{
		x(j) = (b(j) - U(j,j+1)*x(j+1)) / U(j, j);
	}
	return x;
}

// Function 5
// Linear solver using the LU decomposition with no pivoting
// Operation count (2/3) + O(n^2)
VectorXd linear_solve_lu_no_pivoting(const MatrixXd& A, const VectorXd& b)
{
	std::tuple<MatrixXd, MatrixXd> returned_lu;
	returned_lu = lu_no_pivoting(A);
	auto y = forward_subst(std::get<0>(returned_lu), b);
	auto x = backward_subst(std::get<1>(returned_lu), y);
	return x;
}

// Function 6
// Linear solver using LU decomposition with no pivoting for tridiagonal matrices
// Operation count = 8*n + O(1)
VectorXd linear_solve_LU_no_pivoting_tridiag(const MatrixXd& A, const VectorXd& b)
{
	std::tuple<MatrixXd, MatrixXd> returned_lu;
	returned_lu = lu_no_pivoting_tridiag(A);
	auto y = forward_subst_bidiag(std::get<0>(returned_lu), b);
	auto x = backward_subst_bidiag(std::get<1>(returned_lu), y);
	return x;
}

// Function 7
// Linear solver using LU decomposition with row pivoting 
// Operation count for one matrix = (2/3)*n^3 + O(n^2)
// If solving p linear systems, the operation count = (2/3)* p*n^3 + p* O(n^2)
VectorXd linear_solve_lu_row_pivoting(const MatrixXd& A, const VectorXd& b)
{
	std::tuple <MatrixXd, MatrixXd, MatrixXd> returned_lup;
	returned_lup = lu_pivoting(A);
	auto L = std::get<0>(returned_lup);	//The lower triangular matrix is the first element of the tuple
	auto U = std::get<1>(returned_lup); // The upper triangular matrix is the second element of the tuple
	auto P = std::get<2>(returned_lup); // The permutation matrix is the third element of the tuple

	auto y = forward_subst(L, P*b);
	auto x = backward_subst(U, y);
	return x;
}	

// Function 8
// Computing the inverse of a matrix
MatrixXd inverse(const MatrixXd& A)
{
	std::tuple <MatrixXd, MatrixXd, MatrixXd> returned_lup;
	returned_lup = lu_pivoting(A);
	auto L = std::get<0>(returned_lup);	//The lower triangular matrix is the first element of the tuple
	auto U = std::get<1>(returned_lup); // The upper triangular matrix is the second element of the tuple
	auto P = std::get<2>(returned_lup); // The permutation matrix is the third element of the tuple	
	
	int size = A.cols();
	MatrixXd I = MatrixXd::Identity(size, size); // Creating the identity matrix
	MatrixXd inv_matrix(size, size);	// The matrix to be returned
	for (int k = 0; k <= size - 1; ++k)
	{
		auto y = forward_subst(L, P*I.col(k)); // P * k-th column of identity
		auto column_k = backward_subst(U, y);
		inv_matrix.col(k) = column_k;
	}
	return inv_matrix;
}

// Linear solver using cholesky factor
VectorXd linear_solve_cholesky(const MatrixXd& A, const VectorXd& b)
{
    VectorXd x(b.rows()); // x should have same number of rows as b
    MatrixXd choleskyFactor(cholesky(A));

    //First use forward subst to solve Ly = b
    VectorXd y = forward_subst(choleskyFactor.transpose(), b);
    // Now use y to solve Ux = y
    x = backward_subst(choleskyFactor, y);

    return x;
}

// Linear solver using cholesky factor for banded matricies
// This should be switched to use the forward and backward subst of banded matricies
VectorXd linear_solve_banded_cholesky(const MatrixXd& A, const VectorXd& b, int bandSize)
{
	VectorXd x(b.rows()); // x should have same number of rows as b
	MatrixXd choleskyFactor(cholesky_banded_matrix(A, bandSize));

	//First use forward subst to solve Ly = b
	VectorXd y = forward_subst(choleskyFactor.transpose(), b);
	// Now use y to solve Ux = y
	x = backward_subst(choleskyFactor, y);

	return x;
}

#endif // !Solvers_HPP