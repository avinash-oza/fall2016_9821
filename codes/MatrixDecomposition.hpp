#ifndef MatrixDecomposition_HPP
#define MatrixDecomposition_HPP

#include "Eigen/Dense"
#include <tuple>
#include <vector>

using namespace Eigen;

using LU_tuple = std::tuple< MatrixXd, MatrixXd>;
/* LU decomposition without pivoting
		a) Second element of the tuple is the L_no_pivot,
		b) Third element of the tuple is the U_no_pivot,
 */

using LUP_tuple = std::tuple<MatrixXd, MatrixXd, MatrixXd>;
/* LU decomposition with pivoting
		a) Second element of the tuple is the L_no_pivot.
		b) Third element of the tuple is the U_no_pivot.
		c) Fourth element of the tuple is the permutation matrix.
*/

// LU decomposition with no row pivoting
// Operation count = (2/3)* n^3 +O(n^2)
LU_tuple lu_no_pivoting(const MatrixXd& m)
{	// Returns the a tuple where the first element is the L matrix and the second element is the U matrix
	int size = m.cols();

	auto temp = m;	// This way we preserve the elements of A for further use.
	MatrixXd L_no_pivot=  MatrixXd::Zero(size, size);
	MatrixXd U_no_pivot = MatrixXd::Zero(size, size);

	for ( int i = 0; i <= size - 2; ++i)
	{
		for (int k = i; k <= size - 1; ++k)
		{
			U_no_pivot(i,k) = temp(i,k);
			L_no_pivot(k,i) = temp(k,i) / U_no_pivot(i,i);
		}
		for (int j = i + 1; j <= size - 1; ++j)
		{
			for (int k = i + 1; k <= size - 1; ++k)
			{
				temp(j,k) = temp(j,k) - L_no_pivot(j,i) * U_no_pivot(i,k);
			}
		}
	}
	L_no_pivot(size - 1, size - 1) = 1; U_no_pivot(size - 1,size - 1) = temp(size - 1, size - 1);
	LU_tuple myTuple = std::make_tuple(L_no_pivot, U_no_pivot);
	return myTuple;
}

// LU decomposition for trigiagonal matrices NO ROW PIVOTING
// Operation count = 3*n -3
LU_tuple lu_no_pivoting_tridiag(const MatrixXd A)
{
	int size = A.cols();
	auto temp = A;	// This way we preserve the elements of A for further use.
	MatrixXd L_no_pivot = MatrixXd::Zero(size, size);
	MatrixXd U_no_pivot = MatrixXd::Zero(size, size);

	for (int i = 0; i <= size - 1; ++i)
	{
		L_no_pivot(i, i) = 1; 
		L_no_pivot(i + 1, i) = temp(i + 1, i) / temp(i, i);
		U_no_pivot(i, i) = temp(i, i); 
		temp(i + 1, i + 1) = temp(i + 1, i + 1) - L_no_pivot(i + 1, i)*U_no_pivot(i, i + 1);
	}
	LU_tuple myTuple = std::make_tuple(L_no_pivot, U_no_pivot);
	return myTuple;
}

// The lu decomposition with row pivoting
LUP_tuple lu_pivoting(const MatrixXd& m)
{
	int size = m.cols();
	MatrixXi P_vector = MatrixXi::Zero(size,1);			// A matrix of ints of size X
	MatrixXd L_pivot =  MatrixXd::Identity(size, size);	// Initializing the L matrix to Identity matrix
	MatrixXd U_pivot =  MatrixXd::Zero(size, size);		// Initializing the L matrix to Zero matrix
	MatrixXd P_matrix = MatrixXd::Zero(size, size);		// Initializing the permutation matrix to zero

	// Initializing the P_vector to such that it holds the indeces of the permutation matrix.
	for (int i = 0; i < size; i++)
	{
		P_vector(i,0) = i;		// Initializing the P elements to be elements 0:n-1...
	}
	auto temp = m;	

	for (int i = 0; i <= size - 2; i++)
	{
		// We need to find the greatest element in magnitute of the elements of the matrix
		// on the column i, starting from the element with index i until the element with index n-1.

		int index_of_maximum = i; // Assume the maximum is the first element in the range
		for (int k = i + 1; k <= size - 1; k++)
		{
			if (abs(temp(index_of_maximum,i)) < abs(temp(k,i)))
				index_of_maximum = k;
		} // At this point we have the index where the maximum of the elements of column i, starting from the row i to the bottom.

		std::vector<double> vv((size) - i);	// The size of the vector decreases as i increases.
		for (int k = i; k <= size - 1; ++k)
		{
			vv[k - i] = temp(i,k);
			// This is the row that will be switched with the row where the index of the row containing the maximum element.
		}
		for (int k = i; k <= size - 1; ++k)
		{
			temp(i,k) = temp(index_of_maximum,k);
			temp(index_of_maximum,k) = vv[k - i];
		}

		// Altering the rows of the column vector P (which represents the identity matrix).
		int cc = P_vector(i,0);
		P_vector(i,0) = P_vector(index_of_maximum,0);
		P_vector(index_of_maximum,0) = cc;
		if (i > 0)
		{
			// Switching the rows of the L matrix
			std::vector<double> ww(i - 1 + 1);
			for (int k = 0; k <= i - 1; ++k)
			{
				ww[k] = L_pivot(i,k);
				L_pivot(i,k) = L_pivot(index_of_maximum,k);
				L_pivot(index_of_maximum,k) = ww[k];
			}
		}
		for (int j = i; j <= size - 1; ++j)
		{
			U_pivot(i,j) = temp(i,j);					// Compute row i of matrix U
			L_pivot(j,i) = temp(j,i) / U_pivot(i,i);	// Compute column i of matrix L
		}

		for (int j = i + 1; j <= size - 1; ++j)
		{
			for (int k = i + 1; k <= size - 1; ++k)
			{
				temp(j,k) = temp(j,k) - L_pivot(j,i) * U_pivot(i,k);
			}
		}

	}
	L_pivot(size - 1,size - 1) = 1;
	U_pivot(size - 1,size - 1) = temp(size - 1,size - 1);

	// Creating the permutation matrix from the permutation vector
	for (int i = 0; i <= size - 1; ++i)
	{
		P_matrix(i, P_vector(i,0)) = 1;	// Cool... The element of one array is the index of another.
	}
	LUP_tuple myTuple = std::make_tuple(L_pivot, U_pivot, P_matrix);
	return myTuple;
}


// Cholesky decomposition from table 6.1
// Calculates cholesky decomposition. Breaks down if matrix is not SPD
// Operation count is (2/3)*n^3 + O(n^2)
MatrixXd cholesky(const MatrixXd &A)
{
	int originalMatrixSize = A.rows();

    MatrixXd copiedA(A); //copy the matrix A into one that can be changed
    MatrixXd U(originalMatrixSize, originalMatrixSize);

    U.setZero(); // Initialize U to the zero matrix

	int n = A.rows() - 1; // Calculate the location of the last element (0 indexed)
	for (int i = 0; i <= n - 1; i++)
	{
        U(i, i) = std::sqrt(copiedA(i,i));
        for (int k = i + 1; k <= n; k++)
        {
            U(i, k) = copiedA(i, k)/U(i, i); // compute row i of U
        }

        for (int j = i + 1; j <= n; j++)
        {
            for (int k = j; k <= n; k++)
            {
                copiedA(j,k) = copiedA(j,k) - U(i,j)*U(i,k);
            }
        }
	}
    U(n, n) = std::sqrt(copiedA(n,n));

    return U;
}

//Cholesky decomposition for banded matrix
// Assumes the band size is one unless given
MatrixXd cholesky_banded_matrix(const MatrixXd & A, int bandSize)
{
    int originalMatrixSize = A.rows();

    MatrixXd copiedA(A); //copy the matrix A into one that can be changed
    MatrixXd U(originalMatrixSize, originalMatrixSize);

    U.setZero(); // Initialize U to the zero matrix

    int n = A.rows() - 1; // Calculate the location of the last element (0 indexed)
    for (int i = 0; i <= n - 1; i++)
    {
        U(i, i) = std::sqrt(copiedA(i,i));
        for (int k = i + 1; k <= std::min(n, k + bandSize); k++)
        {
            U(i, k) = copiedA(i, k)/U(i, i); // compute row i of U
        }

        for (int j = i + 1; j <= std::min(n, j + bandSize); j++)
        {
            for (int k = j; k <= std::min(n, k + bandSize); k++)
            {
                copiedA(j,k) = copiedA(j,k) - U(i,j)*U(i,k);
            }
        }
    }
    U(n, n) = std::sqrt(copiedA(n,n));

    return U;

}

#endif // !MatrixDecomposition_HPP