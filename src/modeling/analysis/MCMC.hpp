#ifndef RJ_MCMC_HPP
#define RJ_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class Model;
class Move;
class Settings;

/**
 * @brief Markov Chain Monte Carlo using an arbitrary set of moves on a particular Model. 
 */
class MCMC{
    public:
                                MCMC(void)=delete;                                              
                                MCMC(Model* m, std::vector<Move>& mv, Settings& s, bool dBO);   // Constructor
        void                    burnin();                                                       // Run the Markov chain without collecting samples
        void                    run();                                                          // Run the Markov chain and collect the samples to approximate the posterior
    private:
        double                  GibbsIteration(double currentLnPosterior);                      // Run a single iteration of the Markov chain

        BayesianOptimizer       optim;                                                          // The Bayesian optimizer to use if we are optimizing parameters with Bayesian optimization
        bool                    disableBayesOpt;                                                // Are we using Bayesian optimization?
        double                  totalWeight;                                                    // The total weight of all of the MCMC proposals
        int                     bayesOptFreq;                                                   // How often to perform Bayesian optimization
        int                     bayesOptIter;                                                   // How many iterations of Bayesian optimization to perform
        int                     numBurnIn;                                                      // How many burn in iterations to perform
        int                     numIter;                                                        // How many iterations to perform in the sampling run
        int                     printFreq;                                                      // How often to print the Markov chain's state to the screen
        int                     sampleFreq;                                                     // How often to sample the Markov chain
        int                     tuneFreq;                                                       // How often to tune the MCMC proposals
        Model*                  model;                                                          // The model we are running MCMC on
        std::vector<Move>&      moves;                                                          // The proposals we are using during the MCMC
};

#endif