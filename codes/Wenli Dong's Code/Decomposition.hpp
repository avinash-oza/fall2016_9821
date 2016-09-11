#ifndef DECOMPOSITION_HPP
#define DECOMPOSITION_HPP

#include<Eigen/Dense>

void lu_no_pivoting(Eigen::MatrixXd &A);
void lu_row_pivoting(Eigen::MatrixXd &A);

Eigen::MatrixXd cholesky_decomposition(Eigen::MatrixXd &A);

#endif