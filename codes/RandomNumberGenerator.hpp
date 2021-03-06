#ifndef RandomNumberGenerator
#define RandomNumberGenerator
#include <vector>
#include "Eigen/Dense"

using namespace Eigen;

using namespace std;


class UniformVariableGenerationMethod
{
    public:
        virtual VectorXd generateNSamples(int numberOfSamples) = 0;
        virtual void resetGenerator() {};
        virtual double generateSingleSample()  = 0;
};

class LinearCongruentialGenerator : public UniformVariableGenerationMethod
{
    public:
        virtual VectorXd generateNSamples(int numberOfSamples) {
            // Variables needed to avoid overflow
            int q, r, k;
            double a_current;
            double x_current = double(x0);
            VectorXd result=VectorXd::Zero(numberOfSamples);
            for (int i = 0; i < numberOfSamples; ++i)
            {
                result[i] = generateSingleSample();
            }
            return result;
        }

        virtual void resetGenerator()
        {
            // reset the seed so this can be reused
            x0 = 1;
        }

    double generateSingleSample()
    {
        long k = x0 / q;
        x0 = a*(x0-k*q)-k*r;
        if (x0 < 0)
            x0 = x0 + m;
        return (x0+0.0) / m;
    }

    private:
        long x0 = 1;
        const int a = 39373;
        const long int m = pow(2, 31) - 1;
        const int c = 0;
        const long int q = floor(m / a);
        const long int r = m%a;

};

//class AntitheticGenerationMethod: public LinearCongruentialGenerator
//{
//    public:
//        virtual VectorXd generateNSamples(int numberOfSamples) {
//             We take half the samples from U and the other half from 1-U
//            VectorXd uniformSamples = LinearCongruentialGenerator::generateNSamples(numberOfSamples/2);
//            VectorXd difference(uniformSamples.size());
//            difference.setConstant(1.0);
//
//            VectorXd complementUniformSamples = difference - uniformSamples;
//
//            VectorXd returnValues(numberOfSamples);
//            returnValues << uniformSamples, complementUniformSamples;
//            std::cout << complementUniformSamples << std::endl;
//            return returnValues;
//
//        }
//
//};

class NormalVariableGenerationMethod
{
    public:
        virtual VectorXd generateNSamples(int numberOfSamples, UniformVariableGenerationMethod &uniformMethod) = 0;
        double StandardNormalInverse(double u) {
            double a0, a1, a2, a3;
            double b0, b1, b2, b3;
            double c0, c1, c2, c3, c4, c5, c6, c7, c8;
            a0 = 2.50662823884; a1 = -18.61500062529; a2 = 41.39119773534; a3 = 25.44106049637;
            b0 = -8.47351093090;	b1 = 23.08336743743; b2 = -21.06224101826; b3 = 3.13082909833;
            c0 = 0.3374754822726147; c1 = 0.9761690190917186; c2 = 0.1607979714918209;
            c3 = 0.0276438810333863; c4 = 0.0038405729373609; c5 = 0.0003951896511919;
            c6 = 0.0000321767881768; c7 = 0.0000002888167364; c8 = 0.0000003960315187;

            double y,r,x;
            y = u - 0.5;
            if (abs(y) < 0.42){
                r = y*y;
                x = y*(((a3*r + a2)*r + a1)*r + a0) / ((((b3*r + b2)*r + b1)*r + b0)*r + 1);
            }
            else {
                r = u;
                if (y > 0)
                    r = 1 - u;
                r = log(-log(r));
                x = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 + r*c8)))))));
                if (y < 0)
                    x = -x;
            }
            return x;
        }

};


class InverseTransformMethod : public NormalVariableGenerationMethod
{
    public:
        virtual VectorXd generateNSamples(int numberOfSamples, UniformVariableGenerationMethod &uniformMethod) {
            // Create the uniform(0,1) random variable.
            // use linear congruential generator to get samples here
            VectorXd result = uniformMethod.generateNSamples(numberOfSamples);
            for (int i = 0; i < numberOfSamples; ++i)
            {
                result[i] = this->StandardNormalInverse(result[i]);
            }
            return result;
        }
};

class AcceptanceRejectionMethod: public NormalVariableGenerationMethod
{
    public:
    virtual VectorXd generateNSamples(int numberOfSamples, UniformVariableGenerationMethod &uniformMethod) {
        //double u1, u2, u3;

        // If assuming a sample size s, then we need the largest integer that is divisible by 3, that
        // does not exceed size.
        int corrected_size = 3 * (numberOfSamples / 3);

        // If instead we need to generate 3*N sample size, comment out the line above and use the following line:
//        int corrected_size = 3 * numberOfSamples;

        VectorXd UniformSample = uniformMethod.generateNSamples(corrected_size);
        vector<double> std_UniformSample; // To dynamically resize

        for (int i = 0; i <corrected_size; i+=3)
        {
            double x;
            x = -log(UniformSample[i]);
            if (UniformSample[i + 1] > exp(-0.5*(x - 1)*(x - 1))) {}
            else
            {
                if (UniformSample[i + 2] <= 0.5)
                    x = -x;
                std_UniformSample.push_back(x);
            }
        }
        int final_size=std_UniformSample.size();
        VectorXd final_vector = VectorXd::Zero(final_size);
        for (int i = 0; i < final_size; ++i)
            final_vector[i] = std_UniformSample[i];
        return final_vector;
    }
};


//VectorXd LinearCongruentialGenerator(int size){
//
//}




// Acceptance Rejection Method
//VectorXd AcceptanceRejectionMethod(int size) {
//
//}

class BoxMullerMethod: public NormalVariableGenerationMethod
{
    public:
        virtual VectorXd generateNSamples(int numberOfSamples, UniformVariableGenerationMethod &uniformMethod) {
            VectorXd normalDistSamples = VectorXd::Zero(numberOfSamples);

            double u1, u2;
            int i = 0; // how many samples are generated right now
            while (i < numberOfSamples - 1)
            {
                double x, y, z1, z2;
                u1 = 2 * uniformMethod.generateSingleSample() - 1;
                u2 = 2 * uniformMethod.generateSingleSample() - 1;
                x = u1 * u1 + u2 * u2;
                if (x < 1) {
                    y = sqrt(-2 * log(x) / x);
                    z1 = u1 * y;
                    z2 = u2 * y;
                    normalDistSamples[i] = z1;
                    normalDistSamples[i + 1] = z2;
                    i += 2; // only record non rejected samples
                }
            }
            return normalDistSamples;
        }
};


#endif // !RandomNumberGenerator
