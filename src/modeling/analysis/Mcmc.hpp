#ifndef MCMC_HPP
#define MCMC_HPP
#include <vector>

class Model;
class Parameter;
class MoveScheduler;
class EventManager;

class Mcmc{
    public:
        Mcmc(void)=delete;
        Mcmc(Model* m, MoveScheduler* mS);
        void run(int numCycles, int screenIterations, int fileIterations);
    private:
        MoveScheduler* moveScheduler;
        Model* model;
};

#endif