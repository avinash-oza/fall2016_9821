#include <tuple>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <vector>
#include "Eigen/Dense"
#include "IterativeMethods.hpp"
#include "StoppingCriterionSolvers.hpp"
#include "Solvers.hpp"
#include "Examples.hpp"
#include "OptionPricers.hpp"
#include "BlackScholes.hpp"
//#include "Options.hpp"
#include "RandomNumberGenerator.hpp"
#include "MonteCarlo.hpp"
#include "FiniteDifferenceMethods.hpp"


using namespace Eigen;
using namespace std;


//defined at end
MatrixXd generateHW3Matrix(int N);
void runOneIteration(int N, double w);
void printCSVMatrix(std::string stringToPrint, const MatrixXd& myMatrix);
void Question3();
void hw8();

int main() {
    /// Keep this line to make the decimals always print out
    std::cout << std::fixed << std::setprecision(9);
    hw8();
//    Question3();
//    decompositionExamples();
//    verifyCholeskyDecomposition();
//    exam2013();
//    for (int i = 10; i <= 1280; i *= 2 )
//    {
//        calculateTrinomialTreesForN(i);
//    }

    return 0;
}

void hw8()
{
    hw8gLeft  gLeft;
    hw8gRight gRight;
    hw8f f;
    uExact uExact1;

    double xLeft = -2.0;
    double xRight = 2.0;
    double tauFinal = 1.0;
    long M = 2;
    long N = 8;
    cout << setprecision(16);
    double tol = std::pow(10, -6);
    double omega = 1.2;

//    PDESolver solver(gLeft, gRight, f, 0, tauFinal, xLeft, xRight, M, N);
//    MatrixXd fEulerResult = solver.CrankNicolson(LU, tol, omega);
//    MatrixXd bEulerResult = solver.CrankNicolson(LU, tol, omega);
//    printCSVMatrix("Print" ,bEulerResult);

    for (int i = 0; i < 4; ++i)
    {
        M *= 4;
        N *= 2;

        cout << M << "," << N  << ",";
        PDESolver solver(gLeft, gRight, f, 0, tauFinal, xLeft, xRight, M, N);
//        MatrixXd fEulerResult = solver.forwardEuler();
        MatrixXd bEulerResult = solver.backwardEuler(SOR, tol, omega);
//        MatrixXd bEulerResult = solver.CrankNicolson(SOR, tol, omega);
//        std::cout << fEulerResult << std::endl;
        std::cout << solver.MaxPointwiseApproximationError(bEulerResult, uExact1)
        << "," <<  solver.RootMeanSquaredError(bEulerResult, uExact1) << std::endl;
    }


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

void Question3()
{
    // Calculating the option price.
    double spot, strike, interest, vol, maturity, div;
    BlackScholesOption option;
    double price;
    spot = 50; strike = 55; interest = 0.04; vol = 0.3; maturity = 0.5; div = 0;
    price = option.putPrice(spot, strike, maturity, div, interest, 0, vol);

    // Creating the vector containing the value for the number of random variables.
//    int N_vector_size = 10;
    int N_vector_size = 1;
    Eigen::VectorXi N_vector = Eigen::VectorXi::Zero(N_vector_size);
//    N_vector << 1, 2, 4, 8, 16, 32, 64, 128, 256, 512;
    N_vector << 1;
    N_vector *= 10000;


    // Part a: Iverse Transform Method
    InverseTransformMethod inverseTransformMethod;
    LinearCongruentialGenerator uniformMethod;
    MonteCarloMethod monteCarloPricer;

    monteCarloPricer.runMonteCarloForPaths(spot, strike, interest, vol, div, maturity, N_vector, inverseTransformMethod, uniformMethod, price);
//    double vol2 = 0.2, spot2 = 30;
//    spot = 25; strike = 50; interest = 0.05; vol = 0.3; maturity = 0.5; div = 0;
//
//    BasketOptionMonteCarloMethod basketOptionMonteCarloMethod(0.35);
//    basketOptionMonteCarloMethod.runMonteCarloForPaths(spot, spot2, strike, interest, vol, vol2, div, maturity, N_vector, inverseTransformMethod, uniformMethod, price);


//     Part b: Acceptance Rejection Method
//    AcceptanceRejectionMethod acceptanceRejectionMethod;
//    runMonteCarloForPaths(spot, strike, interest, vol, div, maturity, N_vector, acceptanceRejectionMethod, uniformMethod, price);


//    std::cout << std::endl << std::endl;
//     Part c: Box  Muller Method
//    BoxMullerMethod boxMullerMethod;
//    runMonteCarloForPaths(spot, strike, interest, vol, div, maturity, N_vector, boxMullerMethod, uniformMethod, price);
//    std::cout << std::endl << std::endl;
}