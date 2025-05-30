#include "M3S2Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/M3S2Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/M3S2Matrix.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

M3S2Mcmc::M3S2Mcmc(M3S2Model* m, TreeParameter* t, M3S2Matrix* cm, Settings& s) : 
    model(m), codonMatrix(cm), tree(t), generalUpdates(3), stationaryUpdates(5), treeUpdates(0) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;

    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
    tipsLog = s.tipsOutput;
    ancestralLog = s.ancestralStatesOutput;

    kChoice = s.kWeight;
    treeChoice = s.treeWeight + kChoice;
    omega1Choice = s.omegaWeight + treeChoice;
    omega2Choice = s.omegaWeight + omega1Choice;
    omega3Choice = s.omegaWeight + omega2Choice;
    r1Choice = s.rWeight + omega3Choice;
    r2Choice = s.rWeight + r1Choice;
    gammaChoice = s.rWeight + r2Choice;
    stationaryChoice = s.stationaryWeight + gammaChoice;

    treeUpdates = (int)(tree->getTree()->getBranchLengths().size() * 0.5);

    model->regenerateLikelihood();
    model->accept();
}

double M3S2Mcmc::GibbsIteration(double currentLnPosterior){
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
    else if(randomMove < omega1Choice){
        updater = [this]() { return codonMatrix->updateOmega1(); };
    }
    else if(randomMove < omega2Choice){
        updater = [this]() { return codonMatrix->updateOmega2(); };
    }
    else if(randomMove < omega3Choice){
        updater = [this]() { return codonMatrix->updateOmega3(); };
    }
    else if(randomMove < r1Choice){
        updater = [this]() { return codonMatrix->updateR1(); };
    }
    else if(randomMove < r2Choice){
        updater = [this]() { return codonMatrix->updateR2(); };
    }
    else if(randomMove < gammaChoice){
        updater = [this]() { return codonMatrix->updateGamma(); };
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

void M3S2Mcmc::burnin(){
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
            std::cout << "Accept Rates Since Last Tuning Iteration:" << 
                         "\tTree Rate=" << (double)tree->treeAcceptCount/(double)tree->treeCount << 
                         "\tBranch Rate=" << (double)tree->branchAcceptCount/(double)tree->branchCount <<
                         "\tStationary Rate=" << (double)codonMatrix->stationaryAcceptCount/(double)codonMatrix->stationaryCount <<
                         "\tK Rate=" << (double)codonMatrix->kAcceptCount/(double)codonMatrix->kCount <<
                         "\tOmega[1] Rate=" << (double)codonMatrix->omega1AcceptCount/(double)codonMatrix->omega1Count <<
                         "\tOmega[2] Rate=" << (double)codonMatrix->omega2AcceptCount/(double)codonMatrix->omega2Count <<
                         "\tOmega[3] Rate=" << (double)codonMatrix->omega3AcceptCount/(double)codonMatrix->omega3Count <<
                         "\tR[1] Rate=" << (double)codonMatrix->r1AcceptCount/(double)codonMatrix->r1Count <<
                         "\tR[2] Rate=" << (double)codonMatrix->r2AcceptCount/(double)codonMatrix->r2Count <<
                         "\tGamma Rate=" << (double)codonMatrix->gammaAcceptCount/(double)codonMatrix->gammaCount << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }
}

void M3S2Mcmc::run(){
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

    if(tipsLog != ""){
        fs.open(tipsLog, std::ofstream::out);
        fs << model->tipsHeader();
        fs.close();
    }

    if(ancestralLog != ""){
        fs.open(ancestralLog, std::ofstream::out);
        fs << model->ancestralHeader();
        fs.close();
    }

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

            if(tipsLog != "" || ancestralLog != ""){
                auto reconstruction = model->reconstructionOut(n);
    
                if(tipsLog != ""){
                    fs.open(tipsLog, std::ofstream::app);
                    fs << std::get<0>(reconstruction);
                    fs.close();
                    fs.clear();
                }

                if(ancestralLog != ""){
                    fs.open(ancestralLog, std::ofstream::app);
                    fs << std::get<1>(reconstruction);
                    fs.close();
                    fs.clear();
                }
            }
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }
}
