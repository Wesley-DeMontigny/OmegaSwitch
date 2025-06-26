#ifndef SB_MCMC_HPP
#define SB_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class SBModel;
class Parameter;
class TreeParameter;
class SBMatrix;
class Settings;

class SBMcmc{
    public:
        SBMcmc(void)=delete;
        SBMcmc(SBModel* m, TreeParameter* t, SBMatrix* cm, Settings& s, bool dBO);
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
        int truncations;

        double treeChoice;
        double stationaryChoice;
        double kChoice;
        double omegaChoice;
        double rChoice;
        double proportionChoice;

        TreeParameter* tree;
        SBMatrix* codonMatrix;
        SBModel* model;

        BayesianOptimizer optim;
        bool disableBayesOpt;

        std::string analysisLog;
        std::string treeLog;
        std::string tipsLog;
        std::string ancestralLog;
        std::string branchLog;
        
        double GibbsIteration(double currentLnPosterior);
};

#endif