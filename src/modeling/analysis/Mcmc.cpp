#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "MoveScheduler.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <cmath>
#include <iostream>

Mcmc::Mcmc(Model* m, MoveScheduler* mS) : model(m), moveScheduler(mS) {   }

void Mcmc::run(int numCycles, int screenIterations, int fileIterations){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    model->regenerateLikelihood();
    model->accept();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    std::cout << model->tabularHeader() << std::endl;
    for(int n = 1; n <= numCycles; n++){
        double lnProposalRatio = moveScheduler->updateRandom();
        model->regenerateLikelihood();

        double newLnPosterior = model->lnLikelihood() + model->lnPrior();

        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + lnPosteriorRatio;

        if(std::log(rng.uniformRv()) < lnR){
            model->accept();
            currentLnPosterior = newLnPosterior;
        }
        else{
            model->reject();
        }

        if(n % screenIterations == 0){
            std::cout << model->tabularOut(n) << std::endl;
        }
        if(n % fileIterations == 0){

        }
    }
}
