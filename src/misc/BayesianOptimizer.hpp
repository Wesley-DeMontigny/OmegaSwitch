#ifndef BAYESIAN_OPTIMIZER_HPP
#define BAYESIAN_OPTIMIZER_HPP

#include <vector>
#include "Math.hpp"

/**
 * @brief A simple struct for packaging a vector of parameters and
 * an associated score.
 */
struct ParamScorePair {
    std::vector<double> params;
    double score;
    
    // Useful so you can easily rank parameter combinations
    bool operator<(const ParamScorePair& other) const {
        return score < other.score;
    }
};

/**
 * @brief Tunes the parameters used in our MCMC analyses by modeling auto-correlation as a function of the
 * tuning parameters using a Gaussian process with an ARD kernel.
 */
class BayesianOptimizer {
    public:
                                            BayesianOptimizer(void)=delete;
                                            BayesianOptimizer(int nP, int s);                                                               // Constructor

        double                              objective(std::vector<std::vector<double>>& r) {return smoothAverageAutocorrelation(r);}        // Compute the objective that Bayesian optimization seeks to minimize
        std::vector<double>                 getMaximum();                                                                                   // Return the sample corresponding to the highest objective function evalulation
        std::vector<double>                 maximizeAcquisition();                                                                          // Expected improvement acquisition maximization
        std::vector<std::vector<double>>    getMaximumN(int N);                                                                             // Return the top N samples corresponding to the highest objective function evalulation
        void                                registerSample(std::vector<double> s, double o);                                                // Take in sample from MCMC
        void                                setBounds(std::vector<double>& diff, std::vector<double>& mean);                                // Define the bounds of our search using the change of the parameters so far and the current center
        void                                updateGaussianProcess();                                                                        // Update choleskyFactor, alpha, and the length scale hyperparams
    private:
        double                              autocorrelationScore(const std::vector<std::vector<double>>& r, int start, int end);            // Construct the autocorrelation for some subset of the posterior sample
        double                              kernel(std::vector<double>& a, std::vector<double>& b);                                         // ARD kernel
        double                              logPrior();                                                                                     // Compute the log prior for the kernel parameters of the Gaussian process
        double                              marginalLogLikelihood();                                                                        // Compute the marginal log likelihood of the Gaussian process for hyperparameter estimation
        std::vector<double>                 nLLGradient();                                                                                  // Compute the gradient of the marginal log likelihood with respect to the hyperparameters
        std::vector<double>                 nLPGradient();                                                                                  // Compute the gradient of the log prior with respect to the hyperparameters
        double                              smoothAverageAutocorrelation(const std::vector<std::vector<double>>& r);                        // Compute the smoothed autocorrelation score across the whole posterior sample
        double                              UCB(std::vector<double>& sample, int numSamples);                                               // Evaluate the upper confidence bound acquisition function
        void                                updateCholesky();                                                                               // Update the Cholesky factor for the covaraince matrix
        
        double                              sampleMean;                                                                                     // The mean autocorrelation from the MCMC samples
        double                              sampleVariance;                                                                                 // The variance in the autocorrelation from the MCMC samples
        int                                 iterationsPerSample;                                                                            // How many MCMC iterations are associated with each sample
        int                                 numParams;                                                                                      // The size of the optimization space - the number of tunable parameters in our MCMC
        Matrix<double>                      choleskyFactor;                                                                                 // The Cholesky factor for the covariance matrix
        Matrix<double>                      kernelMatrix;                                                                                   // The covariance matrix of the Gaussian process
        std::vector<double>                 alpha;                                                                                          // The alpha vector from  Rasmussen and Williams' Cholesky based marginal log likelihood
        std::vector<double>                 hyperparams;                                                                                    // The length scales for each dimension in the ARD kernel
        std::vector<double>                 lowerBound;                                                                                     // The lower bound for the Latin Hypercube sampler used to initialize the GP update particles
        std::vector<double>                 meanSampleComponents;                                                                           // The mean sample from the MCMC
        std::vector<double>                 upperBound;                                                                                     // The upper bound for the Latin Hypercube sampler use to initialize the GP update particles
        std::vector<ParamScorePair>         samples;                                                                                        // The tunable MCMC parameters associated with a particular batch of MCMC iterations
};

#endif