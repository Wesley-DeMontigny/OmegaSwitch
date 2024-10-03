#ifndef MCMC_HPP
#define MCMC_HPP

class PosteriorNode;
class MoveScheduler;
class EventManager;

class Mcmc{
    public:
        Mcmc(void)=delete;
        Mcmc(PosteriorNode* pN, MoveScheduler* mS);
        void run(int numCycles, EventManager* e);
    private:
        MoveScheduler* moveScheduler;
        PosteriorNode* posterior;
};

#endif