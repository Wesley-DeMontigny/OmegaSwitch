#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/PosteriorNode.hpp"
#include "moves/MoveScheduler.hpp"
#include "moves/Move.hpp"
#include "events/EventManager.hpp"
#include <cmath>
#include <iostream>

Mcmc::Mcmc(PosteriorNode* pN, MoveScheduler* mS) : posterior(pN), moveScheduler(mS) {   }

void Mcmc::run(int numCycles, EventManager* e){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    posterior->regenerate();
    posterior->accept();
    posterior->clean();
    double currentLnPosterior = posterior->lnPosterior();

    for(int n = 1; n <= numCycles; n++){
        Move* m = moveScheduler->getMove();
        double lnProposalRatio = m->update();
        posterior->regenerate();
        double newLnPosterior = posterior->lnPosterior();

        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + lnPosteriorRatio;

        if(std::log(rng.uniformRv()) < lnR){
            m->markAccepted();
            posterior->accept();
            posterior->clean();
            currentLnPosterior = newLnPosterior;
        }
        else{
            m->markRejected();
            posterior->reject();
            posterior->clean();
        }

        e->call(n);
    }
}
