#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

VectorXd forward_subst(const MatrixXd & L, const VectorXd & b)
{
    int lengthOfVector = b.rows();
    VectorXd x(lengthOfVector);

    x(0) = b(0)/(L(0,0));

    for (int j = 1; j < lengthOfVector ; j++)
    {
        double sum = 0;
        // mystery as to why this next line is this
        for (int k = 0; k <= j - 1 ; k++)
        {
            sum += L(j, k) * x(k);
        }
        x(j) = (b(j) - sum)/L(j,j);
    }

    return x;
}

VectorXd backward_subst(const MatrixXd & U, const VectorXd &b)
{
    int lengthOfVector = b.rows();
    VectorXd x(lengthOfVector);

    int n = lengthOfVector - 1;
    x(n) = b(n)/(U(n, n));

    for (int j = n - 1; j >= 0; j--)
    {
        double sum = 0;
        for (int k = j + 1 ; k <= n; k++)
        {
            sum += U(j, k)*x(k);
        }
        x(j) = (b(j) - sum)/U(j,j);
    }
    return x;
}


int main()
{
//    MatrixXd L(4,4);
//    VectorXd b(4);
//
//    L << 100, 0,0,0,
//            6, 106, 0, 0,
//            8,8, 108, 0,
//            5,5,5,105;
//    b << 98, 104, 111,102;
//
//    std::cout << L << std::endl;
//    std::cout << b << std::endl;
//    std::cout << forward_subst(L, b) << std::endl
    MatrixXd U(3,3);
    VectorXd b(3);

    U << 1,2,3,0,4,5,0,0,7;
    b << 1,2,3;

    std::cout << U << std::endl;
    std::cout << b << std::endl;
    std::cout << backward_subst(U, b) << std::endl;
}