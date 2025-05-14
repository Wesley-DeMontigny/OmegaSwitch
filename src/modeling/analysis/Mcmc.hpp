#ifndef MCMC_HPP
#define MCMC_HPP
#include <vector>
#include <string>

class Model;
class Parameter;
class TreeParameter;
class CodonMultiMatrix;
class DirichletProcessPrior;
class Settings;

class Mcmc{
    public:
        Mcmc(void)=delete;
        Mcmc(Model* m, TreeParameter* t, CodonMultiMatrix* cmm, DirichletProcessPrior* dpp, Settings& s);
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
        CodonMultiMatrix* codonMatrix;
        DirichletProcessPrior* dpp;

        std::string analysisLog;
        std::string treeLog;
        std::string dppLog;
        std::string tipsLog;

        Model* model;
        
        double GibbsIteration(double currentLnPosterior);
};

#endif