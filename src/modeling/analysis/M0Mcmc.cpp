#include "M0Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/M0Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/M0Matrix.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

M0Mcmc::M0Mcmc(M0Model* m, TreeParameter* t, M0Matrix* cm, Settings& s) : 
    model(m), codonMatrix(cm), tree(t), generalUpdates(3), stationaryUpdates(5), treeUpdates(0) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;

    kChoice = s.kWeight;
    treeChoice = s.treeWeight + kChoice;
    omegaChoice = s.omegaWeight + treeChoice;
    stationaryChoice = s.stationaryWeight + omegaChoice;

    treeUpdates = (int)(tree->getTree()->getBranchLengths().size() * 0.5);

    model->regenerateLikelihood();
    model->accept();
}

double M0Mcmc::GibbsIteration(double currentLnPosterior){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double randomMove = rng.uniformRv() * stationaryChoice;

    std::function<double()> updater;
    int numUpdates = generalUpdates;

    if(randomMove < kChoice){
        updater = [this]() { return codonMatrix->updateK(); };
    }
    else if(randomMove < treeChoice){
        updater = [this]() { return tree->update(); };
        numUpdates = treeUpdates;
    }
    else if(randomMove < omegaChoice){
        updater = [this]() { return codonMatrix->updateOmega(); };
    }
    else if(randomMove < stationaryChoice){
        updater = [this]() { return codonMatrix->updateStationary(); };
        numUpdates = stationaryUpdates;
    }

    for(int i = 0; i < numUpdates; i++){
        double lnProposalRatio = updater();
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

    return currentLnPosterior;
}

void M0Mcmc::burnin(){
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
            std::cout << "Accept Rates Since Last Tuning Iteration:" << 
                         "\tTree Rate=" << (double)tree->treeAcceptCount/(double)tree->treeCount << 
                         "\tBranch Rate=" << (double)tree->branchAcceptCount/(double)tree->branchCount <<
                         "\tStationary Rate=" << (double)codonMatrix->stationaryAcceptCount/(double)codonMatrix->stationaryCount <<
                         "\tK Rate=" << (double)codonMatrix->kAcceptCount/(double)codonMatrix->kCount <<
                         "\tOmega Rate=" << (double)codonMatrix->omegaAcceptCount/(double)codonMatrix->omegaCount << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }
}

void M0Mcmc::run(){
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
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }
}
