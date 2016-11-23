#include <tuple>
#include <iomanip>
#include<iostream>
#include<fstream>
#include <cmath>
#include <time.h>
#include <vector>
#include "Eigen/Dense"
#include "IterativeMethods.hpp"
#include "StoppingCriterionSolvers.hpp"
#include "Solvers.hpp"
#include "OptionPricers.hpp"
#include "BlackScholes.hpp"
//#include "Options.hpp"
#include "FiniteDifferenceMethods.hpp"


using namespace Eigen;
using namespace std;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ",", "\n");
//defined at end
/*
MatrixXd generateHW3Matrix(int N);
void runOneIteration(int N, double w);
void printCSVMatrix(std::string stringToPrint, const MatrixXd& myMatrix);
void Question3();
*/
void hw8();
int main() {
    /// Keep this line to make the decimals always print out
    std::cout << std::fixed << std::setprecision(9);
	/*
	std::ofstream myfile("output1.csv");
	myfile<<(U).format(CSVFormat)<<std::endl;
	myfile.close();
	*/
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
    //long M = 2;
    //long N = 8;
	long M = 4;
    cout << setprecision(16);
    double tol = std::pow(10, -6);
    double omega = 1.2;

	double sigma = 0.35;
	double S0 = 41;
	double q = 0.02;
	double K = 40;
	double T = 0.75;
	double r = 0.04;
	double alphatemp = 0.45;

    hw8fOption fOption(sigma, S0, q, K, T, r);
    hw8gLeftOption gLeftOption(sigma, S0, q, K, T, r);
    hw8gRightOption gRightOption;

//    uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0, double tFinal,
//            double S0, double K, double T, double q, double r, double sigma, int M, double alphatemp) :
//    PDESolver(gLeftFunc, gRightFunc, f, 0, T, 0, 0, 0, 0), S0(S0)
//    EuropeanPutPDESolver(uFunction &gLeftFunc, uFunction &gRightFunc, uFunction &f, double t0,
//            double S0, double K, double T, double q, double r, double sigma, int M, double alphatemp) :
    EuropeanPutPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
//    std::cout << solver.N << std::endl;
//    std::cout << solver.M << std::endl;
    std::cout << solver.CrankNicolson(SOR, tol, omega) << std::endl;
//	double a = (r-q)/(sigma*sigma) - 0.5;
//	double b = pow((r-q)/(sigma*sigma) + 0.5, 2) + 2*q/(sigma*sigma);
//	double xLeftOption = log(S0/K) + (r-q-sigma*sigma/2) * T - 3*sigma*sqrt(T);
//	double xRightOption = log(S0/K) + (r-q-sigma*sigma/2) * T + 3*sigma*sqrt(T);
//	double tauFinalOption = T*sigma*sigma/2;

	
//    PDESolver solver(gLeft, gRight, f, 0, tauFinal, xLeft, xRight, M, N);
//    MatrixXd fEulerResult = solver.CrankNicolson(LU, tol, omega);
//    MatrixXd bEulerResult = solver.CrankNicolson(LU, tol, omega);
//    printCSVMatrix("Print" ,bEulerResult);
	/*
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
	*/
	BlackScholesOption option(S0, K, T, q, r, sigma);
	double Vexact = option.putPrice();
	//PDESolver solver(gLeftOption, gRightOption, fOption, 0, tauFinalOption, xLeftOption, xRightOption, M, alphatemp);
	//MatrixXd EulerResultOption = solver.forwardEuler();
	//MatrixXd EulerResultOption = solver.backwardEuler(LU, tol, omega);
	//MatrixXd EulerResultOption = solver.backwardEuler(SOR, tol, omega);
	//MatrixXd EulerResultOption = solver.CrankNicolson(LU, tol, omega);
	//MatrixXd EulerResultOption = solver.CrankNicolson(SOR, tol, omega);
	//cout<<fEulerResultOption;
	/*
	std::ofstream myfile("output1.csv");
	myfile<<(EulerResultOption).format(CSVFormat)<<std::endl;
	myfile.close();
	*/
	MatrixXd result=MatrixXd::Zero(4,9);
	for(int count=0; count<4; ++count) {
        M *= 4;
//        EuropeanPutPDESolver solver(gLeftOption, gRightOption, fOption, 0, S0, K, T, q, r, sigma, M, alphatemp);
//        std::cout << solver.CrankNicolson(LU, tol, omega) << std::endl;
    }
		// Pointwise Convergence
//		PDESolver solver(gLeftOption, gRightOption, fOption, 0, tauFinalOption, xLeftOption, xRightOption, M, alphatemp);
		//MatrixXd EulerResultOption = solver.forwardEuler();
		//MatrixXd EulerResultOption = solver.backwardEuler(LU, tol, omega);
		//cout<<EulerResultOption;
		//MatrixXd EulerResultOption = solver.backwardEuler(SOR, tol, omega);
		//MatrixXd EulerResultOption = solver.CrankNicolson(LU, tol, omega);
//		MatrixXd EulerResultOption = solver.CrankNicolson(SOR, tol, omega);
//		cout<<count<<endl<<endl;
//		long N = solver.getN();
//		double x=log(S0/K);
//		int i=0;
//		for(i;i<N+1;++i)
//		{
//			if(x<solver.mesh.getX(i))
//				break;
//		}
//		i=i-1;
//		double S1 = K*exp(solver.mesh.getX(i));
//		double S2 = K*exp(solver.mesh.getX(i+1));
//		double V1 = EulerResultOption(M,i)*exp(-a*solver.mesh.getX(i)-b*tauFinalOption);
//		double V2 = EulerResultOption(M,i+1)*exp(-a*solver.mesh.getX(i+1)-b*tauFinalOption);
//		double Uappro = ((solver.mesh.getX(i+1)-x)*EulerResultOption(M,i)+(x-solver.mesh.getX(i))*EulerResultOption(M,i+1))/(solver.mesh.getX(i+1)-solver.mesh.getX(i));
//		double Vappro1 = ((S2-S0)*V1 + (S0-S1)*V2)/(S2-S1);
//		double Vappro2 = exp(-a*x-b*tauFinalOption)*Uappro;


//		double error1=abs(Vappro1-Vexact);
//		double error2=abs(Vappro2-Vexact);
//		result(count,0)=error1;
//		result(count,2)=error2;
		//cout<<error1<<' '<<error2<<endl<<endl;

		//RMS Error
//		VectorXd Vappro = VectorXd::Zero(N+1);
//		VectorXd S = VectorXd::Zero(N+1);
//		double Nrms = 0;
//		double sum = 0;
//		for(int j = 0; j < N+1; ++j)
//		{
//			S(j) = K*exp(solver.mesh.getX(j));
//			Vappro(j) = EulerResultOption(M,j)*exp(-a*solver.mesh.getX(j)-b*tauFinalOption);
//            option.setS(S(j));
//			double price = option.putPrice();
//			if(price > 0.00001*S0)
//			{
//				++Nrms;
//				sum = sum + pow(Vappro(j) - price , 2)/pow(price,2);
//			}
//		}
//		result(count,4)=sqrt(sum/Nrms);
//		if(count>0)
//		{
//			result(count,1) = result(count, 0)/result(count-1,0);
//			result(count,3) = result(count, 2)/result(count-1,2);
//			result(count,5) = result(count, 4)/result(count-1,4);
//		}
//
//		 Greeks
//		double delta = (Vappro(i+1) - Vappro(i))/(S(i+1) - S(i));
//		result(count,6)=delta;
//		cout<<delta;
//		double delta2 = (Vappro(i+2) - Vappro(i+1))/(S(i+2) - S(i+1));
//		double delta0 = (Vappro(i) - Vappro(i-1))/(S(i) - S(i-1));
//		double gamma = (delta2-delta0)/((S(i+2)+S(i+1)-S(i)-S(i-1))/2);
//		result(count,7)=gamma;
//		cout<<gamma;
//		double dT = 2*(solver.mesh.getT(M) - solver.mesh.getT(M-1))/(sigma*sigma);
//		double V1dT = EulerResultOption(M-1,i)*exp(-a*solver.mesh.getX(i)-b*solver.mesh.getT(M-1));
//		double V2dT = EulerResultOption(M-1,i+1)*exp(-a*solver.mesh.getX(i+1)-b*solver.mesh.getT(M-1));
//		double Vt = ((S2-S0)*V1dT + (S0-S1)*V2dT)/(S2-S1);
//		double theta = (Vappro1 - Vt)/dT;
//		result(count,8)=theta;
//		cout<<theta;
//	}
//	cout<<'1';
//	std::ofstream myfile("output2.csv");
//	myfile<<(result).format(CSVFormat)<<std::endl;
//	myfile.close();
//
}


