#ifndef MonteCarlo_HPP
#define MonteCarlo_HPP

#include "RandomNumberGenerator.hpp"
#include <cmath>
#include <string>
#include <tuple>
#include <iostream>


using namespace std;
tuple<double,long int> MonteCarlo(double Spot, double Strike, double Interest, double Volatility,
                                  double Dividend, double Maturity, int NumberOfPaths, NormalVariableGenerationMethod &normalVariableGenerationMethod, UniformVariableGenerationMethod &uniformMethod);

class MonteCarloMethod
{
    public:
        virtual double calculatePriceFromSimulatedPaths(const VectorXd &spotPrices, const VectorXd &disc_payoff_vector, double discountFactor, double S0)
        {
            return disc_payoff_vector.mean();
        }

        virtual void runMonteCarloForPaths(double Spot, double Strike, double Interest, double Volatility, double Dividend, double Maturity,
                                       VectorXi pathsToRunFor, NormalVariableGenerationMethod &transformMethod,
                                       UniformVariableGenerationMethod &uniformMethod, double price)
            {
                for (int i = 0; i < pathsToRunFor.size(); ++i)
                {
                    std::cout << setprecision(12) << setw(5);
                    tuple<double, long int> MonteCarloTuple = MonteCarlo(Spot, Strike, Interest, Volatility, Dividend, Maturity, pathsToRunFor[i], transformMethod, uniformMethod);
                    double monte_carlo_price = std::get<0>(MonteCarloTuple); // Returns the price
                    long int number_simulations = std::get<1>(MonteCarloTuple); // Returns the number of simulations
                    cout << pathsToRunFor[i] << "\t" << number_simulations << "\t" << monte_carlo_price << "\t" << abs(price - monte_carlo_price) << std::endl;
                }
                std::cout << std::endl << std::endl;
            }

        virtual VectorXd adjustSpotPrice(const VectorXd &spotPrices, double discountFactor, double spotPrice)
        {
            return spotPrices;
        }

        virtual tuple<double,long int> MonteCarlo(double Spot, double Strike, double Interest, double Volatility,
                                          double Dividend, double Maturity, int NumberOfPaths, NormalVariableGenerationMethod &normalVariableGenerationMethod, UniformVariableGenerationMethod &uniformMethod)
        {
            VectorXd sample_random_var = normalVariableGenerationMethod.generateNSamples(NumberOfPaths, uniformMethod);
            double discountFactor = exp(-Interest*Maturity);

            // At this point, we have the vector with the sample variables
            // We first find its size.
            int size = sample_random_var.size();

            // Create the vector with the spot prices corresponding with each element in the random sample
            VectorXd spot_price_vector = VectorXd::Zero(size);
            for (int i = 0; i < size; ++i)
            {
                spot_price_vector[i] = calculateLogNormalSpotPrice(Spot, Interest, Volatility, Dividend, Maturity, sample_random_var[i]);
            }

            spot_price_vector = adjustSpotPrice(spot_price_vector, discountFactor, Spot);

            // Create the vector with the discounted payoffs
            VectorXd disc_payoff_vector = VectorXd::Zero(size);
            for (int i = 0; i < size; ++i)
            {
                disc_payoff_vector[i] = exp(-Interest*Maturity)*max(Strike - spot_price_vector[i], 0.0);
            }
            // Then, return the mean of the payoffs
            double optionPrice = calculatePriceFromSimulatedPaths(spot_price_vector, disc_payoff_vector, discountFactor, Spot);
            long int sample_size = size;
            tuple<double, long int> results = make_tuple(optionPrice, sample_size);
            return results;
        }

    virtual double calculateLogNormalSpotPrice(double Spot, double Interest, double Volatility, double Dividend, double Maturity, double normalVariableValue) {
        return Spot * exp((Interest - Dividend - pow(Volatility, 2) / 2.0) * Maturity + Volatility * sqrt(Maturity) * normalVariableValue); }
};

class ControlVariateMonteCarloMethod: public MonteCarloMethod
{
    public:
        virtual double
            calculatePriceFromSimulatedPaths(const VectorXd &spotPrices, const VectorXd &disc_payoff_vector, double discountFactor, double S0) {

                VectorXd S_hat(spotPrices.size());
                S_hat.setConstant(spotPrices.mean());

                VectorXd V_hat(spotPrices.size());
                V_hat.setConstant(disc_payoff_vector.mean());

                VectorXd dfVector(spotPrices.size());
                dfVector.setConstant(1/discountFactor*S0);

                VectorXd meanAdjustedSpotPrices = spotPrices - S_hat;
                VectorXd meanAdjustedPayoffs = disc_payoff_vector - V_hat;

                double b_hat = meanAdjustedSpotPrices.dot(meanAdjustedPayoffs)/(meanAdjustedSpotPrices.dot(meanAdjustedSpotPrices));

                VectorXd W = disc_payoff_vector - b_hat*(spotPrices - dfVector);
                return W.mean();

            }
};

class MomentMatchingMonteCarloMethod: public MonteCarloMethod
{
    public:
            virtual VectorXd adjustSpotPrice(const VectorXd &simulatedPrices, double discountFactor, double spotPrice) {
                double S_hat = simulatedPrices.mean();
                double factorToAdjust = spotPrice/(discountFactor*S_hat);
                return simulatedPrices*factorToAdjust;
        }
};

class MomentMatchingAndControlVariateMonteCarloMethod : public MomentMatchingMonteCarloMethod, ControlVariateMonteCarloMethod
{
    public:
        virtual void runMonteCarloForPaths(double Spot, double Strike, double Interest, double Volatility, double Dividend,
                                           double Maturity, VectorXi pathsToRunFor,
                                           NormalVariableGenerationMethod &transformMethod,
                                           UniformVariableGenerationMethod &uniformMethod, double price) override {
            MomentMatchingMonteCarloMethod::runMonteCarloForPaths(Spot, Strike, Interest, Volatility, Dividend, Maturity, pathsToRunFor,
                                                    transformMethod, uniformMethod, price);
        }

        virtual double calculatePriceFromSimulatedPaths(const VectorXd &spotPrices, const VectorXd &disc_payoff_vector,
                                                        double discountFactor, double S0) override {
            return ControlVariateMonteCarloMethod::calculatePriceFromSimulatedPaths(spotPrices, disc_payoff_vector, discountFactor, S0);
        }
};


class BasketOptionMonteCarloMethod : public MonteCarloMethod
{
    public:
        BasketOptionMonteCarloMethod(double rho): rho(rho) {};

    virtual void runMonteCarloForPaths(double Spot, double Spot2, double Strike, double Interest, double Volatility, double Volatility2, double Dividend, double Maturity,
                                       VectorXi pathsToRunFor, NormalVariableGenerationMethod &transformMethod,
                                       UniformVariableGenerationMethod &uniformMethod, double price)
    {
        for (int i = 0; i < pathsToRunFor.size(); ++i)
        {
            std::cout << setprecision(12) << setw(5);
            tuple<double, long int> MonteCarloTuple = MonteCarlo(Spot, Spot2, Strike, Interest, Volatility, Volatility2, Dividend, Maturity, pathsToRunFor[i], transformMethod, uniformMethod);
            double monte_carlo_price = std::get<0>(MonteCarloTuple); // Returns the price
            long int number_simulations = std::get<1>(MonteCarloTuple); // Returns the number of simulations
            cout << pathsToRunFor[i] << "\t" << number_simulations << "\t" << monte_carlo_price << "\t" << abs(price - monte_carlo_price) << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    virtual tuple<double,long int> MonteCarlo(double Spot, double Spot2, double Strike, double Interest, double Volatility, double Volatility2,
                                          double Dividend, double Maturity, int NumberOfPaths, NormalVariableGenerationMethod &normalVariableGenerationMethod, UniformVariableGenerationMethod &uniformMethod)
        {
            VectorXd sample_random_var = normalVariableGenerationMethod.generateNSamples(2*NumberOfPaths, uniformMethod);
            double discountFactor = exp(-Interest*Maturity);

            // At this point, we have the vector with the sample variables
            // We first find its size.
            int size = sample_random_var.size();
            int vectorSize = NumberOfPaths / 2; // The size of the spot, stock and option price vectors

            // Create the vector with the spot prices corresponding with each element in the random sample
            VectorXd spot_price_vector = VectorXd::Zero(vectorSize);
            VectorXd second_spot_price_vector = VectorXd::Zero(vectorSize);
            for (int i = 0; i < vectorSize; ++i)
            {
                double &firstSample = sample_random_var[2 * i];
                double &secondSample = sample_random_var[2 * i + 1];
                spot_price_vector[i] = calculateLogNormalSpotPrice(Spot, Interest, Volatility, Dividend, Maturity, firstSample);
                // calculate out the final z that should go in for pricing
                double finalStandardNormValue = rho * firstSample + sqrt(1 - rho * rho) * secondSample;
                second_spot_price_vector[i] = calculateLogNormalSpotPrice(Spot2, Interest, Volatility2, Dividend, Maturity, finalStandardNormValue);
            }

//            spot_price_vector = adjustSpotPrice(spot_price_vector, discountFactor, Spot);
            VectorXd strikeVector = VectorXd::Ones(vectorSize);
            strikeVector *= Strike; // Make the vector of strikes
            VectorXd undiscountedOptionValues = spot_price_vector + second_spot_price_vector - strikeVector;


            // Create the vector with the discounted payoffs
            VectorXd disc_payoff_vector = VectorXd::Zero(vectorSize);
            for (int i = 0; i < vectorSize; ++i)
            {
                disc_payoff_vector[i] = discountFactor*max(undiscountedOptionValues[i], 0.0);
            }
            // Then, return the mean of the payoffs
            double optionPrice = calculatePriceFromSimulatedPaths(spot_price_vector, disc_payoff_vector, discountFactor, Spot);
            long int sample_size = size;
            tuple<double, long int> results = make_tuple(optionPrice, sample_size);
            return results;
        }
    private:
        double rho;
};


#endif // !MonteCarlo_HPP
