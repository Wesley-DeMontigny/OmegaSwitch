#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#include <string>
#include <vector>

/**
 * @brief The struct containing all of the settings a user can set from the command line
 * to run our software
 */
struct Settings {
    
                            Settings(int argc,  char* argv[]);
                            Settings(void) = delete;

    void                    print();
    void                    usage();

    bool                    M0 = false;                     // Do inference under a standard codon phylogenetic mode
    bool                    CMM = false;                     // Do inference under a Markov-modulated model with simulataneous inference for the number of evolutionary regimes
    bool                    DPCMM = false;                  // Do inference under an infinite mixture of Markov-modulated models with simulateous inference for the number of evolutionary regimes
    bool                    sequentialTuningSim = false;    // Do a simulation that first does Bayesian optimization for tunable parameters and then does classic tuning
    bool                    simulateM0 = false;             // Simulate under a standard codon model
    bool                    simulateCMM = false;             // Simulate under a Markov-modulated model with 1-5 evolutionary regimes
    bool                    simulateDPCMM = false;          // Simulate under a mixture of Markov-modulated models with 1-3 evolutionary regimes
    double                  dppWeight = 2.0;                // Weight associated with performing Neal's algorithm 8
    double                  expectedCat = 2.0;              // The expected number of categories under the Chinese restaurant process for CMM-DPP
    double                  mcmcmcBeta = 1.0;               // Inverse temperature for the tempered chain; values below 1.0 enable MCMCMC
    double                  kLambda = 3.0;                  // The rate parameter for the transition/transversion exponential prior
    double                  kWeight = 1.0;                  // Weight associated with updating the transition/transversion rate parameter
    double                  omegaLambda = 3.0;              // The rate parameter for the non-synonymous/synonymous exponential prior
    double                  omegaWeight = 1.0;              // Weight associated with updating the non-synonymous/synonymous rate parameter
    double                  rjWeight = 1.0;                 // Weight associated with the Reversible-Jump (Metropolis-Hastings-Green) move
    double                  rLambda = 3.0;                  // The rate parameter for the evolutionary regime swapping exponential prior
    double                  rWeight = 1.0;                  // Weight associated with updating the rate of swapping evolutionary rates
    double                  stationaryWeight = 2.0;         // Weight associated with updating the stationary distribution of the rate matrix
    double                  treeLengthMean = 35.0;          // The mean parameter associated with the tree length prior
    double                  treeLengthSD = 2.5;             // The SD parameter associated with the tree length prior
    double                  treeWeight = 2.0;               // Weight associated with updating the tree topology or branch lengths
    int                     bayesOpt = 0;                   // How many iterations of Bayesian optimization to perform on the tunable parameters
    int                     bayesOptFrequency = 0;          // How many MCMC iterations a single iteration of Bayesian optimization uses
    int                     burnInIterations = 1000;        // How many MCMC iterations to discard before beginning sampling
    int                     numGibbs = 25;                  // How many sites to sample in a single iteration of Neal's algorithm 8
    int                     numIterations = 5000;           // How many MCMC iterations to sample for
    int                     numSimulations = 0;             // How many simulations to perform under the specified model
    int                     mcmcmcSwapFrequency = 10;      // How many iterations between attempted swaps in MCMCMC
    int                     printFrequency = 10;            // How many MCMC iterations between printing the program state to the screen
    int                     sampleFrequency = 10;           // How many MCMC iterations between sampling events
    id_t                    threads = 15;                    // How many threads to use in the analysis
    int                     tuneFrequency = 1000;            // How many MCMC iterations between tuning events during the burn-in

    std::string             ancestralStatesOutput = "";     // The file to output all ancestral states to
    std::string             branchOutput = "";              // The file to output the branch length trace to
    std::string             dppOutput = "";                 // The file to output the DPP trace to
    std::string             tree = "";                      // The tree topology to enforce during MCMC
    std::string             mcmcOutput = "";                // The file to output the general analysis trace to
    std::string             nexusInput = "";                // The nexus input file
    std::string             simulationOutput = "";          // The file to output the true simulated states to
    std::string             tipsOutput = "";                // The file to output the tip dN/dS trace to
    std::string             treeOutput = "";                // The file to output the newick trace to
};

#endif
