#ifndef RJ_MCMC_HPP
#define RJ_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class RJModel;
class Parameter;
class TreeParameter;
class RJMatrix;
class Settings;

class RJMcmc{
    public:
        RJMcmc(void)=delete;
        RJMcmc(RJModel* m, TreeParameter* t, RJMatrix* cm, Settings& s, bool dBO);
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
        double rChoice;
        double rjChoice;

        TreeParameter* tree;
        RJMatrix* codonMatrix;
        RJModel* model;

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