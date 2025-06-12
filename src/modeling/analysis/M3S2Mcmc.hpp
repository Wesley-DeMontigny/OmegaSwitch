#ifndef M3S2_MCMC_HPP
#define M3S2_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class M3S2Model;
class Parameter;
class TreeParameter;
class M3S2Matrix;
class Settings;

class M3S2Mcmc{
    public:
        M3S2Mcmc(void)=delete;
        M3S2Mcmc(M3S2Model* m, TreeParameter* t, M3S2Matrix* cm, Settings& s);
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

        TreeParameter* tree;
        M3S2Matrix* codonMatrix;
        M3S2Model* model;
        BayesianOptimizer optim;

        std::string analysisLog;
        std::string treeLog;
        std::string tipsLog;
        std::string ancestralLog;
        std::string branchLog;
        
        double GibbsIteration(double currentLnPosterior);
};

#endif