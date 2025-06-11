#ifndef BAYESIAN_OPTIMIZER_HPP
#define BAYESIAN_OPTIMIZER_HPP

#include <vector>
#include "Math.hpp"

class BayesianOptimizer {
    public:
        BayesianOptimizer(void)=delete;
        BayesianOptimizer(int nP, int s);
        ~BayesianOptimizer();

        static double objective(std::vector<double> r);
        std::vector<double> maximizeAcquisition(); // TODO: Implement expected improvement acquisition function
        std::vector<double> getMaximum(); // Return the sample corresponding to the highest objective function evalulation
        void updateGaussianProcess(); // Update choleskyFactor, alpha, and the length scale hyperparams

        void registerSample(std::vector<double> s, double o) {samples.push_back(s); objectives.push_back(o);}
        void setBounds(std::vector<double>& diff, std::vector<double>& mean, std::vector<bool>& sign); // Define the bounds of our search using the change of the parameters so far and the current center
    private:
        int numParams;
        int iterationsPerSample;

        std::vector<std::vector<double>> samples; // The MCMC parameters associated with a particular batch of MCMC iterations
        std::vector<double> objectives; // The values of the objective function for each of the samples

        std::vector<double> hyperparams; // The length scales of the ARD kernel
        std::vector<double> upperBound;
        std::vector<double> lowerBound;
        std::vector<double> initialSampleScale;  

        Matrix<double> choleskyFactor;
        std::vector<double> alpha;

        static double autocorrelationScore(const std::vector<double>& r, int end); //Mahendran et al. 2010
        double kernel(std::vector<double>& a, std::vector<double>& b); // ARD kernel
        double marginalLogLikelihood();
        void updateCholesky();
};

#endif