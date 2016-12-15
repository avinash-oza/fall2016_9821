//
// Created by Tushar Chawla on 12/9/16.
//

#include "Exams.hpp"


int main()
{
    std::cout.precision(9);
    Exam2015::Question2::Question2Part3();
    //std::cout << std::endl << std::endl;
    //Exam2015::Question2::Question2Part4();
    //std::cout << std::endl << std::endl;
    //Exam2015::Question2::Question2Part5();
    //std::cout << std::endl << std::endl;
    //Exam2015::Question3::Question3part1();
   // Exam2015::Question1::Question1Part3();

    //Exam2015::Question3::Question3Part3();

   /* BlackScholesPutOption blackScholesPutOption(47.00287094, 50.0, (7.0/12.0)-0.003054101,  0.01, 0.03, 0.25);
    std::cout << std::endl;
    std::cout << blackScholesPutOption.price() << std::endl;*/

   /* double S = 41.0;
    double K = 40.0;
    double T = 1.0; // make sure to do double division
    double q = 0.01;
    double r = 0.03;
    double sigma = 0.25; //needed for cnstruction, not used however
    OPTION_TYPE optionType = PUT;

    std::cout.precision(6);
    AmericanBinomialTreePricer americanBinomialTreePricer(S,K,T,q,r,sigma,optionType);
    //TREE_RESULT tree_result = americanBinomialTreePricer.calculateTree(190);
    std::cout << "Final Implied Vol " << std::endl;
    std::cout<< americanBinomialTreePricer.calculateImpliedVolatilityViaSecantMethod(190,2.0,.1,.7,0.000001) << std::endl;
    */
}
