#ifndef MonteCarlo_HPP
#define MonteCarlo_HPP

#include "RandomNumberGenerator.hpp"
#include <cmath>
#include <string>
#include <tuple>
#include <iostream>


using namespace std;
tuple<double,long int> MonteCarlo(double Spot, double Strike, double Interest, double Volatility, double Dividend, double Maturity, int NumberOfPaths, RandomNumberGenerationMethod &transformMethod);


void runMonteCarloForPaths(double Spot, double Strike, double Interest, double Volatility, double Dividend, double Maturity, VectorXi pathsToRunFor, RandomNumberGenerationMethod &transformMethod, double price)
{
    for (int i = 0; i < pathsToRunFor.size(); ++i)
    {
        std::cout << setprecision(12) << setw(5);
        tuple<double, long int> MonteCarloTuple = MonteCarlo(Spot, Strike, Interest, Volatility, Dividend, Maturity, pathsToRunFor[i], transformMethod);
        double monte_carlo_price = std::get<0>(MonteCarloTuple); // Returns the price
        long int number_simulations = std::get<1>(MonteCarloTuple); // Returns the number of simulations
        cout << pathsToRunFor[i] << "\t" << number_simulations << "\t" << monte_carlo_price << "\t" << abs(price - monte_carlo_price) << std::endl;
    }
    std::cout << std::endl << std::endl;
}

tuple<double,long int> MonteCarlo(double Spot, double Strike, double Interest, double Volatility, double Dividend, double Maturity, int NumberOfPaths, RandomNumberGenerationMethod &transformMethod)
{
    VectorXd sample_random_var = transformMethod.generateNSamples(NumberOfPaths);

    // At this point, we have the vector with the sample variables
	// We first find its size.
	int size = sample_random_var.size();

	// Create the vector with the spot prices corresponding with each element in the random sample
	VectorXd spot_price_vector = VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		spot_price_vector[i] = Spot*exp((Interest - Dividend - pow(Volatility, 2)/2.0)*Maturity + Volatility*sqrt(Maturity)*sample_random_var[i]);
	}

	// Create the vector with the discounted payoffs
	VectorXd disc_payoff_vector = VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		disc_payoff_vector[i] = exp(-Interest*Maturity)*max(Strike - spot_price_vector[i], 0.0);
	}
	// Then, return the mean of the payoffs
	double average = disc_payoff_vector.mean();
	long int sample_size = size;
	tuple<double, long int> results = make_tuple(average, sample_size);
	return results;
}
#endif // !MonteCarlo_HPP
