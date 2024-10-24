#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "MoveScheduler.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

Mcmc::Mcmc(Model* m, MoveScheduler* mS) : model(m), moveScheduler(mS) {   }

void Mcmc::burnin(int numCycles, int screenIterations, int tuningIterations){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    model->regenerateLikelihood();
    model->accept();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numCycles; n++){
        if(n % screenIterations == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
        }
        if(n % tuningIterations == 0){
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

void Mcmc::run(int numCycles, int screenIterations, int fileIterations){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    model->regenerateLikelihood();
    model->accept();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    std::string tabularHeader = model->tabularHeader();
    std::cout << tabularHeader;

    std::string analysisLog = "C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/analysis.log";
    std::string treeLog = "C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/trees.log";
    std::string dppLog = "C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/dpp.log";
    std::fstream fs;
    fs.open(analysisLog, std::fstream::out);
    fs << tabularHeader;
    fs.close();
    fs.clear();
    fs.open(treeLog, std::fstream::out);
    fs << model->treeHeader();
    fs.close();
    fs.clear();
    fs.open(dppLog, std::fstream::out);
    fs << model->dppHeader();
    fs.close();
    fs.clear();

    for(int n = 1; n <= numCycles; n++){
        if(n % screenIterations == 0){
            std::cout << model->tabularOut(n);
        }
        if(n % fileIterations == 0){
            fs.open(analysisLog, std::fstream::app);
            fs << model->tabularOut(n);
            fs.close();
            fs.clear();
            fs.open(treeLog, std::fstream::app);
            fs << model->treeOut(n);
            fs.close();
            fs.clear();
            fs.open(dppLog, std::fstream::app);
            fs << model->dppOut(n);
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
