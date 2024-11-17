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
    tipsLog = s.tipsOutput;
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

    std::ofstream fs;
    fs.open(analysisLog, std::ofstream::out);
    fs << tabularHeader;
    fs.close();

    fs.open(treeLog, std::ofstream::out);
    fs << model->treeHeader();
    fs.close();

    fs.open(dppLog, std::ofstream::out);
    fs << model->dppHeader();
    fs.close();

    fs.open(tipsLog, std::ofstream::out);
    fs << model->tipsHeader();
    fs.close();

    for(int n = 1; n <= numIter; n++){
        if(n % printFreq == 0){
            std::cout << model->tabularOut(n);
        }
        if(n % sampleFreq == 0){
            fs.open(analysisLog, std::ofstream::app);
            fs << model->tabularOut(n);
            fs.close();
            fs.clear();

            fs.open(treeLog, std::ofstream::app);
            fs << model->treeOut(n);
            fs.close();
            fs.clear();

            fs.open(dppLog, std::ofstream::app);
            fs << model->dppOut(n);
            fs.close();
            fs.clear();

            model->reconstructTips();
            fs.open(tipsLog, std::ofstream::app);
            fs << model->tipsOut(n);
            fs.close();
            fs.clear();
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
