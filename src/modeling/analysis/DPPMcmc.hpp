#ifndef DPP_MCMC_HPP
#define DPP_MCMC_HPP
#include <vector>
#include <string>

class DPPModel;
class Parameter;
class TreeParameter;
class DPPMatrix;
class DirichletProcessPrior;
class Settings;

class DPPMcmc{
    public:
        DPPMcmc(void)=delete;
        DPPMcmc(DPPModel* m, TreeParameter* t, DPPMatrix* cm, DirichletProcessPrior* dpp, Settings& s);
        void burnin();
        void run();
    private:
        int numIter;
        int numBurnIn;
        int printFreq;
        int tuneFreq;
        int sampleFreq;
        int generalUpdates;
        int stationaryUpdates;
        int treeUpdates;

        double treeChoice;
        double stationaryChoice;
        double kChoice;
        double rChoice;
        double dppChoice;
        double omegaChoice;

        TreeParameter* tree;
        DPPMatrix* codonMatrix;
        DirichletProcessPrior* dpp;
        DPPModel* model;

        std::string analysisLog;
        std::string treeLog;
        std::string dppLog;
        std::string tipsLog;
        std::string ancestralLog;
        
        double GibbsIteration(double currentLnPosterior);
};

#endif