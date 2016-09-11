#include"Decomposition.hpp"
#include"Solvers.hpp"
#include<Eigen/Dense>
#include<iostream>

void lu_no_pivoting(Eigen::MatrixXd &A)
{
	// Initialize L with an identity matrix and U with a square matrixs with all the entries equal to 0
	Eigen::MatrixXd L=Eigen::MatrixXd::Identity(A.rows(),A.cols());
	Eigen::MatrixXd U=Eigen::ArrayXXd::Zero(A.rows(),A.cols());

	for(int j=0;j<A.rows()-1;++j)
	{
		for(int k=j;k<A.rows();++k)
		{
			U(j,k)=A(j,k); // Calculate the jth row of U
			L(k,j)=A(k,j)/U(j,j); // Calculate the jth column of L
		}
		A.block(j+1,j+1,A.rows()-1-j,A.rows()-1-j)=A.block(j+1,j+1,A.rows()-1-j,A.rows()-1-j)
			-L.block(j+1,j,A.rows()-1-j,1)*U.block(j,j+1,1,A.rows()-1-j);
		// Update matrix A[j:n-1,j:n-1]
	}
	U(A.rows()-1,A.cols()-1)=A(A.rows()-1,A.cols()-1); // Calculate U[n-1,n-1]

	std::cout<<L<<std::endl<<std::endl;
	std::cout<<U<<std::endl<<std::endl;

}

void lu_row_pivoting(Eigen::MatrixXd &A)
{
	// Initialize L and P with an identity matrix and U with a square matrixs with all the entries equal to 0
	Eigen::MatrixXd L=Eigen::MatrixXd::Identity(A.rows(),A.cols());
	Eigen::MatrixXd U=Eigen::ArrayXXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd P=Eigen::MatrixXd::Identity(A.rows(),A.cols());

	for(int j=0;j<A.rows()-1;++j)
	{
		// find the row index of the maximal absolute value within A[j:n-1,j] 
		int j_max=j;
		for(int k=j+1;k<A.rows();++k)
		{
			if(abs(A(j_max,j))<abs(A(k,j)))
				j_max=k;
		}

		// Exchange j_max row and the jth row with column from j to n-1 of A
		Eigen::MatrixXd temp=A.block(j,j,1,A.rows()-j);
		A.block(j,j,1,A.rows()-j)=A.block(j_max,j,1,A.rows()-j);
		A.block(j_max,j,1,A.rows()-j)=temp;

		// Exchange j_max row and the jth row of P
		Eigen::MatrixXd temp1=P.block(j,0,1,A.rows());
		P.block(j,0,1,A.rows())=P.block(j_max,0,1,A.rows());
		P.block(j_max,0,1,A.rows())=temp1;

		// Exchange the solved part of j_max row and j_th row of L
		if(j>0)
		{
			Eigen::MatrixXd temp3=L.block(j,0,1,j);
			L.block(j,0,1,j)=L.block(j_max,0,1,j);
			L.block(j_max,0,1,j)=temp3;
		}

		for(int k=j;k<A.rows();++k)
		{
			U(j,k)=A(j,k); // Calculate the jth row of U
			L(k,j)=A(k,j)/U(j,j); // Calculate the jth column of L
		}
		A.block(j+1,j+1,A.rows()-1-j,A.rows()-1-j)=A.block(j+1,j+1,A.rows()-1-j,A.rows()-1-j)
			-L.block(j+1,j,A.rows()-1-j,1)*U.block(j,j+1,1,A.rows()-1-j);
		// Update matrix A[j:n-1,j:n-1]
	}
	U(A.rows()-1,A.rows()-1)=A(A.rows()-1,A.rows()-1); // Calculate U[n-1,n-1]

	std::cout<<P<<std::endl<<std::endl;
	std::cout<<L<<std::endl<<std::endl;
	std::cout<<U<<std::endl<<std::endl;

}

Eigen::MatrixXd cholesky_decomposition(Eigen::MatrixXd &A)
{
	// Initialize U with a square matrixs with all the entries equal to 0
	Eigen::MatrixXd U=Eigen::ArrayXXd::Zero(A.rows(),A.cols());
	for(int j=0;j<A.rows()-1;++j)
	{
		U(j,j)=sqrt(A(j,j));
		for(int k=j+1;k<A.cols();++k)
		{
			U(j,k)=A(j,k)/U(j,j);
		}
		for(int i=j+1;i<A.rows();i++)
			for(int k=i;k<A.cols();k++)
			{
				A(i,k)=A(i,k)-U(j,i)*U(j,k);
			}
	}
	U(A.rows()-1,A.cols()-1)=sqrt(A(A.rows()-1,A.cols()-1));
	return U;


}