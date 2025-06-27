#ifndef BAYESIAN_OPTIMIZER_HPP
#define BAYESIAN_OPTIMIZER_HPP

#include <vector>
#include "Math.hpp"

struct ParamScorePair {
    std::vector<double> params;
    double score;

    bool operator<(const ParamScorePair& other) const {
        return score < other.score;
    }
};

class BayesianOptimizer {
    public:
        BayesianOptimizer(void)=delete;
        BayesianOptimizer(int nP, int s);
        ~BayesianOptimizer();

        double objective(std::vector<std::vector<double>>& r) {return smoothAverageAutocorrelation(r);}
        std::vector<double> maximizeAcquisition(); // Expected improvement acquisition maximization
        std::vector<double> getMaximum(); // Return the sample corresponding to the highest objective function evalulation
        std::vector<std::vector<double>> getMaximumN(int N); // Return the top N samples corresponding to the highest objective function evalulation
        void updateGaussianProcess(); // Update choleskyFactor, alpha, and the length scale hyperparams

        void registerSample(std::vector<double> s, double o);
        void setBounds(std::vector<double>& diff, std::vector<double>& mean); // Define the bounds of our search using the change of the parameters so far and the current center
    private:
        int numParams;
        int iterationsPerSample;

        std::vector<ParamScorePair> samples; // The MCMC parameters associated with a particular batch of MCMC iterations
        double sampleVariance;

        std::vector<double> hyperparams; // The length scales
        std::vector<double> upperBound;
        std::vector<double> lowerBound;

        Matrix<double> choleskyFactor;
        Matrix<double> kernelMatrix;
        std::vector<double> alpha;

        double kernel(std::vector<double>& a, std::vector<double>& b); // Linear kernel
        double marginalLogLikelihood();
        std::vector<double> nLLGradient();
        double autocorrelationScore(const std::vector<std::vector<double>>& r, int start, int end);
        double smoothAverageAutocorrelation(const std::vector<std::vector<double>>& r);
        double averageProportionalJumpingDistance(const std::vector<std::vector<double>>& r);
        double expectedImprovement(std::vector<double>& sample, double currentMaxObjective, int numSamples);
        double UCB(std::vector<double>& sample, int numSamples);
        void updateCholesky();
};

#endif