#ifndef M0_MCMC_HPP
#define M0_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class M0Model;
class Parameter;
class TreeParameter;
class M0Matrix;
class Settings;

class M0Mcmc{
    public:
        M0Mcmc(void)=delete;
        M0Mcmc(M0Model* m, TreeParameter* t, M0Matrix* cm, Settings& s, bool dBO);
        void burnin();
        void run();
    private:
        int numIter;
        int numBurnIn;
        int printFreq;
        int tuneFreq;
        int sampleFreq;
        int bayesOptFreq;
        int bayesOptIter;
        int generalUpdates;
        int stationaryUpdates;
        int treeUpdates;

        double treeChoice;
        double stationaryChoice;
        double kChoice;
        double omegaChoice;

        TreeParameter* tree;
        M0Matrix* codonMatrix;
        M0Model* model;

        BayesianOptimizer optim;
        bool disableBayesOpt;

        std::string analysisLog;
        std::string treeLog;
        std::string branchLog;
        
        double GibbsIteration(double currentLnPosterior);
};

#endif