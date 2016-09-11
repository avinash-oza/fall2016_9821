#include "Solvers.hpp"
#include"Decomposition.hpp"
#include<Eigen/Dense>
#include<iostream>
#include<fstream>

void main()
{
	// Matrix multiplication
	/*
	Eigen::MatrixXd A(2,10);
	A<<7.91,8.01,-3.77,-9.43,7.63,-9.25,-8.98,-4.3,9.87,-9.07,
		-2.97,-3.39,7.47,3.85,-8.9,-1.24,-3.23,9.25,7.23,-2.5;
	   
	Eigen::MatrixXd B(10,1);
	B<<-0.73,-6.01,9.76,-0.05,3.29,-2.49,-9.16,4.3,-9.78,-0.96;
	   
	   
	Eigen::MatrixXd M=matrix_multiply(A,B);
	std::cout<<M<<std::endl;
	*/
	
	
	// Forward substitutions 
	/*
	Eigen::MatrixXd L(7,7);
	L<<-1.6,0,0,0,0,0,0,
       4.86,3.83,0,0,0,0,0,
	   -0.28,1.91,1.93,0,0,0,0,
	   -2.43,2.39,4.55,3.4,0,0,0,
	   2.64,-1.83,-0.37,2.79,1.28,0,0,
	   -4.03,2.12,2.74,-0.68,0.92,0.38,0,
	   4.24,-2.05,3.73,-0.42,-1.75,0.65,2.81;
	Eigen::VectorXd b(7);
	b<<49.274,491.2186,-371.223,-373.504,-243.644,-214.96,281.6024;
	Eigen::VectorXd x = forward_subst(L,b);
	std::cout<<x<<std::endl;
	*/
    	

	//Solve a linear system using forward substitution
	/*
	Eigen::MatrixXd L(3,3);
	L<<1,0,0,
       4,5,0,
	   7,8,9;
	Eigen::MatrixXd b(3,2);
	b<<0,3,
	   1,4,
	   2,5;
	Eigen::MatrixXd x= linearsystem_forward_subst(L,b);
	std::cout<<x<<std::endl;
	*/

	//Backward substitutions
	/*
	Eigen::MatrixXd U(5,5);
	U<<0.03,0.74,-0.75,0.63,0.45,
       0,0.54,0.86,0.66,0.53,
	   0,0,-0.99,-0.27,-0.84,
	   0,0,0,0.59,-0.32,
	   0,0,0,0,-0.99;
	Eigen::VectorXd b(5);
	b<<16.4315,15.42426,-4.32748,13.31876,5.36982;
	Eigen::VectorXd x= backward_subst(U,b);
	std::cout<<x<<std::endl;
	*/

	//Solve a linear system using backward substitution
	/*
	Eigen::MatrixXd U(4,4);
	U<<3,-1,2,-1,
       0,2,-1,3,
	   0,0,4,2,
	   0,0,0,1;
	Eigen::MatrixXd b(4,2);
	b<<0.333333,1.66667,
	   1.166667,3.8333,
	   0.875,1.875,
	   -0.916667,-5.58333;
	Eigen::MatrixXd x= linearsystem_backward_subst(U,b);
	std::cout<<x<<std::endl;
	*/
	
	
	//LU decomposition without pivoting
	/*
	Eigen::MatrixXd A(6,6);
	A<<0.88, -0.04, -0.24, 0.69, -0.83,0.97,
	  0.23, -0.19, -0.23, 0.73, -0.88,-0.32,
	  0.23,-0.22, -0.65, -0.73, 0.84,0.08,
	  -0.21,0.56,-0.94,-0.97,0.31,-0.65,
	  0.08,-0.44,-0.77, -0.68, 0.98,-0.54,
	  0.91,-0.27,0.16,0.98, 0.76,-0.67;
	lu_no_pivoting(A);
	*/

	
	//LU decomposition with row pivoting
	/*
	Eigen::MatrixXd A(6,6);
	A<<0.62, -0.04,-0.24,0.69, -0.83,0.97,
	   0.23, -0.19, -0.23, 0.73, -0.88,-0.32,
	  0.23,-0.22, -0.65,  -0.73,0.84,0.08,
	   -0.21,0.56, -0.94,-0.97,  0.31,-0.65,
	   0.08, -0.44,-0.77, -0.68, 0.08,-0.54,
	   0.91,-0.27,0.16,0.98,0.76,-0.67;
	lu_row_pivoting(A);
	*/

	
	// Cholesky decomposition
	/*
	Eigen::MatrixXd A(4,4);
	A<< 9,-3, 6,-3,
	   -3, 5,-4, 7,
	    6,-4,21, 3,
		-3, 7, 3,15;
	Eigen::MatrixXd U=cholesky_decomposition(A);
	std::cout<<U<<std::endl;
	*/


    //Solve a linear system using choleskey decomposition    
    /*
	Eigen::MatrixXd A(4,4);
	A<< 9,-3, 6,-3,
	   -3, 5,-4, 7,
	    6,-4,21, 3,
		-3, 7, 3,15;
	Eigen::MatrixXd b(4,2);
	b<<1,5,
	   2,6,
	   3,7,
	   4,8;
	Eigen::MatrixXd x=linearsystem_choleskey_decomposition(A,b);
	std::cout<<x<<std::endl;
	*/

   // Jocobi iteration, Guass Siedel iteration, SOR iteration
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

   //Eigen::VectorXd x=Jacobi_Iteration(A,b,x_0,tol);

   //Eigen::VectorXd x=GaussSiedel_Iteration(A,b,x_0,tol);

   /*
   double w=1.15;
   Eigen::VectorXd x= SOR_Iteration(A,b,x_0,tol,w);
   */

   std::ofstream myfile; 
   myfile.open("iv_Output.txt");
   for(double w=1.02; w<2.00; w=w+0.02)
   {
	   myfile<<SOR_Iteration(A,b,x_0,tol,w)<<std::endl;
	   //std::cout<<b-A*x<<std::endl<<std::endl;
   }
   myfile.close();
   
   //std::cout<<x<<std::endl<<std::endl;
   //std::cout<<b-A*x<<std::endl;
   
   
   

 

    
}