#ifndef MCMC_HPP
#define MCMC_HPP
#include <vector>
#include <string>

class Model;
class Parameter;
class MoveScheduler;
class EventManager;
class Settings;

class Mcmc{
    public:
        Mcmc(void)=delete;
        Mcmc(Model* m, MoveScheduler* mS, Settings& s);
        void burnin();
        void run();
    private:
        int numIter;
        int numBurnIn;
        int printFreq;
        int tuneFreq;
        int sampleFreq;
        std::string analysisLog;
        std::string treeLog;
        std::string dppLog;
        MoveScheduler* moveScheduler;
        Model* model;
};

#endif