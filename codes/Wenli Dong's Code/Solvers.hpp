#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include<Eigen/Dense>



// Direct methods
Eigen::VectorXd backward_subst(const Eigen::MatrixXd &U, const Eigen::VectorXd &b);
Eigen::VectorXd forward_subst(const Eigen::MatrixXd &L, const Eigen::VectorXd &b);

Eigen::MatrixXd linearsystem_backward_subst(const Eigen::MatrixXd &U, const Eigen::MatrixXd &b);
Eigen::MatrixXd linearsystem_forward_subst(const Eigen::MatrixXd &U, const Eigen::MatrixXd &b);

Eigen::MatrixXd matrix_multiply(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);

Eigen::MatrixXd linearsystem_choleskey_decomposition(Eigen::MatrixXd &A, const Eigen::MatrixXd &b);

// Iterative methods
double norm_2(const Eigen::VectorXd &x);
Eigen::VectorXd Jacobi_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x_0, double tol);
Eigen::VectorXd GaussSiedel_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::VectorXd &x_0, double tol);
int SOR_Iteration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::VectorXd &x_0, double tol, double w);




#endif