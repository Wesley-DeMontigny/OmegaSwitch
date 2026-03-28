#ifndef RJ_MCMC_HPP
#define RJ_MCMC_HPP
#include <vector>
#include <string>
#include "misc/BayesianOptimizer.hpp"

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
                                MCMC(Model* coldModel, std::vector<Move>& coldMoves,
                                     Model* temperedModel, std::vector<Move>& temperedMoves,
                                     Settings& s, bool dBO);                                     // Constructor for two-chain MCMCMC
        void                    burnin();                                                       // Run the Markov chain without collecting samples
        void                    run();                                                          // Run the Markov chain and collect the samples to approximate the posterior
    private:
        double                  GibbsIteration(Model* activeModel, std::vector<Move>& activeMoves,
                                               double currentLnPosterior, double beta,
                                               double currentTotalWeight);                      // Run a single iteration of the Markov chain
        double                  getBetaForChain(int chainIndex) const;                         // Return the inverse temperature for a chain
        Model*                  getSampleModel() const;                                         // Return the chain currently targeting the posterior
        Model*                  getTemperedModel() const;                                       // Return the tempered chain
        double                  initializeChain(int chainIndex);                                // Initialize a chain and return its raw log-posterior
        void                    printBurninAcceptanceRates() const;                             // Print move acceptance diagnostics
        bool                    attemptSwap();                                                  // Attempt a temperature swap between the two chains

        BayesianOptimizer       optim;                                                          // The Bayesian optimizer to use if we are optimizing parameters with Bayesian optimization
        bool                    disableBayesOpt;                                                // Are we using Bayesian optimization?
        bool                    hasTemperedChain;                                               // Whether we are running two-chain MCMCMC
        double                  hotBeta;                                                        // The inverse temperature of the tempered chain
        int                     bayesOptFreq;                                                   // How often to perform Bayesian optimization
        int                     bayesOptIter;                                                   // How many iterations of Bayesian optimization to perform
        int                     numBurnIn;                                                      // How many burn in iterations to perform
        int                     numIter;                                                        // How many iterations to perform in the sampling run
        int                     printFreq;                                                      // How often to print the Markov chain's state to the screen
        int                     sampleFreq;                                                     // How often to sample the Markov chain
        int                     sampleChainIndex;                                               // The chain currently assigned to the posterior
        int                     swapFreq;                                                       // How often to attempt chain swaps
        int                     swapAcceptCount;                                                // Number of accepted swaps
        int                     swapCount;                                                      // Number of proposed swaps
        int                     tuneFreq;                                                       // How often to tune the MCMC proposals
        double                  totalWeights[2];                                                // The total move weights for each chain
        double                  currentLnPosteriors[2];                                         // The current raw log-posteriors of the chains
        Model*                  models[2];                                                      // The models we are running MCMC on
        std::vector<Move>*      moves[2];                                                       // The proposals we are using during the MCMC
};

#endif
