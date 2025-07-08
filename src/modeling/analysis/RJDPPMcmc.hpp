#ifndef RJ_DPP_MCMC_HPP
#define RJ_DPP_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class RJDPPModel;
class Parameter;
class TreeParameter;
class RJDPPMatrix;
class RJDirichletProcessPrior;
class Settings;

class RJDPPMcmc{
    public:
        RJDPPMcmc(void)=delete;
        RJDPPMcmc(RJDPPModel* m, TreeParameter* t, RJDPPMatrix* cm, RJDirichletProcessPrior* dpp, Settings& s, bool dBO);
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
        double rChoice;
        double dppChoice;
        double omegaChoice;
        double rjChoice;

        TreeParameter* tree;
        RJDPPMatrix* codonMatrix;
        RJDirichletProcessPrior* dpp;
        RJDPPModel* model;

        BayesianOptimizer optim;
        bool disableBayesOpt;

        std::string analysisLog;
        std::string treeLog;
        std::string dppLog;
        std::string tipsLog;
        std::string ancestralLog;
        std::string branchLog;
        
        double GibbsIteration(double currentLnPosterior);
};

#endif