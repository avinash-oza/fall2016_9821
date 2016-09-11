#include"Decomposition.hpp"
#include"Solvers.hpp"
#include<Eigen/Dense>
#include<iostream>

// Direct methods
Eigen::VectorXd backward_subst(const Eigen::MatrixXd &U, const Eigen::VectorXd &b)
{
	Eigen::VectorXd x(U.rows());
	x(U.rows()-1)=b(U.rows()-1)/U(U.rows()-1,U.rows()-1);
	for(int j=U.rows()-2;j>=0;--j)
	{
		double sum=0;
		for(int k=j+1;k<U.rows();++k)
		{
			sum=sum+U(j,k)*x(k);
		}
		x(j)=(b(j)-sum)/U(j,j);
	}
	return x;
}

Eigen::VectorXd forward_subst(const Eigen::MatrixXd &L, const Eigen::VectorXd &b)
{

	Eigen::VectorXd x(L.rows());
	x(0)=b(0)/L(0,0);
	for(int j=1;j<L.rows();++j)
	{
		double sum=0;
		for(int k=0;k<j;++k)
		{
			sum=sum+L(j,k)*x(k);
		}
		x(j)=(b(j)-sum)/L(j,j);
	}
	return x;
}

Eigen::MatrixXd linearsystem_backward_subst(const Eigen::MatrixXd &U, const Eigen::MatrixXd &b)
{
	Eigen::MatrixXd X=Eigen::ArrayXXd::Zero(b.rows(),b.cols());
	for(int k=0;k<X.cols();++k)
	{
		X.block(0,k,X.rows(),1)=backward_subst(U,b.block(0,k,b.rows(),1));
	}
	return X;

}

Eigen::MatrixXd linearsystem_forward_subst(const Eigen::MatrixXd &L, const Eigen::MatrixXd &b)
{
	Eigen::MatrixXd X=Eigen::ArrayXXd::Zero(b.rows(),b.cols());
	for(int k=0;k<X.cols();++k)
	{
		X.block(0,k,X.rows(),1)=forward_subst(L,b.block(0,k,b.rows(),1));
	}
	return X;

}

Eigen::MatrixXd matrix_multiply(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
{
	Eigen::MatrixXd M=Eigen::ArrayXXd::Zero(A.rows(),B.cols());
	for(int j=0;j<A.rows();++j)
		for(int k=0;k<B.cols();++k)
		{
			// Inner product function: v1.dot(v2)
			Eigen::VectorXd v1=A.block(j,0,1,A.cols()).transpose();
			Eigen::VectorXd v2=B.block(0,k,B.rows(),1);
			M(j,k)=v1.dot(v2);
		}
	return M;

}

Eigen::MatrixXd linearsystem_choleskey_decomposition(Eigen::MatrixXd &A, const Eigen::MatrixXd &b)
{
	Eigen::MatrixXd U=cholesky_decomposition(A);
	const Eigen::MatrixXd Y= linearsystem_forward_subst(U.transpose(),b);
	Eigen::MatrixXd X=linearsystem_backward_subst(U,Y);
	return X;
}

// Iterative methods
double norm_2(const Eigen::VectorXd &x)
{
	return sqrt(x.dot(x));
}

Eigen::VectorXd Jacobi_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x_0, double tol)
{
	Eigen::MatrixXd L=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd U=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd D_Inverse=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	for(int j=0; j<A.rows(); ++j)
		for(int k=0; k<A.cols(); ++k)
		{
			if (j==k)
				D_Inverse(j,k)=1/A(j,k);
			else if (j<k)
				U(j,k)=A(j,k);
			else
				L(j,k)=A(j,k);

		}
	
		Eigen::VectorXd x=x_0;


		Eigen::VectorXd b_new=D_Inverse*b;
		int ic=0;

		// Consecutive approximation criterion
		
	    x_0=Eigen::VectorXd::Zero(A.rows());
		while(norm_2(x-x_0)>tol)
		{
			x_0=x;
			x=-D_Inverse*(L*x+U*x)+b_new;
			ic=ic+1;
		}
		

		//Residue-Based criterion
		/*
		Eigen::VectorXd r_0=b-A*x_0;
		Eigen::VectorXd r=r_0;
		double stop_iter_resid=tol*norm_2(r_0);
		while(norm_2(r)>stop_iter_resid)
		{
			x=-D_Inverse*(L*x+U*x)+b_new;
			r=b-A*x;
			ic=ic+1;
		}
		*/

		std::cout<<ic<<std::endl<<std::endl;
		//std::cout<<x<<std::endl<<std::endl;
		return x;		
}

Eigen::VectorXd GaussSiedel_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x_0, double tol)
{
	Eigen::MatrixXd L=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd U=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd D=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	for(int j=0; j<A.rows(); ++j)
		for(int k=0; k<A.cols(); ++k)
		{
			if (j==k)
				D(j,k)=A(j,k);
			else if (j<k)
				U(j,k)=A(j,k);
			else
				L(j,k)=A(j,k);

		}
	
		int ic=0;
		Eigen::VectorXd x=x_0;
		Eigen::VectorXd b_new=forward_subst(L+D,b);
		Eigen::MatrixXd R=linearsystem_forward_subst(L+D,U);

		// Consecutive approximation criterion
		/*
	    x_0=Eigen::VectorXd::Zero(A.rows());
		while(norm_2(x-x_0)>tol)
		{
		    x_0=x;
			x=-R*x+b_new;
			ic=ic+1;
		}
		*/

		//Residue-Based criterion
		
		Eigen::VectorXd r_0=b-A*x_0;
		Eigen::VectorXd r=r_0;
		double stop_iter_resid=tol*norm_2(r_0);
		while(norm_2(r)>stop_iter_resid)
		{
			x=-R*x+b_new;
			r=b-A*x;
			ic=ic+1;
		}
		

		std::cout<<ic<<std::endl<<std::endl;
		//std::cout<<x<<std::endl<<std::endl;
		return x;		
}

int SOR_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::VectorXd &x_0, double tol, double w)
{
	Eigen::MatrixXd L=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd U=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd D_Inverse=Eigen::MatrixXd::Zero(A.rows(),A.cols());
	Eigen::MatrixXd I=Eigen::MatrixXd::Identity(A.rows(),A.cols());

	for(int j=0; j<A.rows(); ++j)
		for(int k=0; k<A.cols(); ++k)
		{
			if (j==k)
				D_Inverse(j,k)=1/A(j,k);
			else if (j<k)
				U(j,k)=A(j,k);
			else
				L(j,k)=A(j,k);

		}
	
		int ic=0;
		Eigen::VectorXd x=x_0;
		Eigen::VectorXd b_new=w*forward_subst(I+w*D_Inverse*L,D_Inverse*b);
		Eigen::MatrixXd R=linearsystem_forward_subst(I+w*D_Inverse*L,(1-w)*I-w*D_Inverse*U);

		// Consecutive approximation criterion
		
		/*
	    Eigen::VectorXd x0=Eigen::VectorXd::Zero(A.rows());
		while(norm_2(x-x0)>tol)
		{
		    x0=x;
			x=R*x+b_new;
			ic=ic+1;
		}
		*/

		//Residue-Based criterion
		
		Eigen::VectorXd r_0=b-A*x_0;
		Eigen::VectorXd r=r_0;
		double stop_iter_resid=tol*norm_2(r_0);
		while(norm_2(r)>stop_iter_resid)
		{
			x=R*x+b_new;
			r=b-A*x;
			ic=ic+1;
		}
		

		//std::cout<<ic<<std::endl<<std::endl;
		return ic;	
		//return x
}