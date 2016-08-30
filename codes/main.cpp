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


int main()
{
    MatrixXd L(4,4);
    VectorXd b(4);

    L << 100, 0,0,0,
            6, 106, 0, 0,
            8,8, 108, 0,
            5,5,5,105;
    b << 98, 104, 111,102;

    std::cout << L << std::endl;
    std::cout << b << std::endl;
    std::cout << forward_subst(L, b) << std::endl;
}