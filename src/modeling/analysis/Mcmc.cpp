#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "MoveScheduler.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

Mcmc::Mcmc(Model* m, MoveScheduler* mS, Settings& s) : model(m), moveScheduler(mS) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
    dppLog = s.dppOutput;
}

void Mcmc::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    model->regenerateLikelihood();
    model->accept();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

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
    }
}

void Mcmc::run(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    model->regenerateLikelihood();
    model->accept();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    std::string tabularHeader = model->tabularHeader();
    std::cout << tabularHeader;

    std::fstream fsAnalysis;
    fsAnalysis.open(analysisLog, std::fstream::out);
    fsAnalysis << tabularHeader;
    fsAnalysis.close();

    std::fstream fsTree;
    fsTree.open(treeLog, std::fstream::out);
    fsTree << model->treeHeader();
    fsTree.close();

    std::fstream fsDPP;
    fsDPP.open(dppLog, std::fstream::out);
    fsDPP << model->dppHeader();
    fsDPP.close();

    for(int n = 1; n <= numIter; n++){
        if(n % printFreq == 0){
            std::cout << model->tabularOut(n);
        }
        if(n % sampleFreq == 0){
            fsAnalysis.open(analysisLog, std::fstream::app);
            fsAnalysis << model->tabularOut(n);
            fsAnalysis.close();

            fsTree.open(treeLog, std::fstream::app);
            fsTree << model->treeOut(n);
            fsTree.close();

            fsDPP.open(dppLog, std::fstream::app);
            fsDPP << model->dppOut(n);
            fsDPP.close();
        }

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
    }

    fsAnalysis.close();
    fsDPP.close();
    fsTree.close();
}
