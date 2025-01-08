#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/CodonMultiMatrix.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

Mcmc::Mcmc(Model* m, TreeParameter* t, CodonMultiMatrix* cmm, DirichletProcessPrior* d, Settings& s) : 
    model(m), dpp(d), codonMatrix(cmm), tree(t), metropolisUpdates(10) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
    dppLog = s.dppOutput;
    tipsLog = s.tipsOutput;

    kChoice = s.kWeight;
    rChoice = s.rWeight + kChoice;
    omegaChoice = s.omegaWeight + rChoice;
    treeChoice = s.treeWeight + omegaChoice;
    stationaryChoice = s.stationaryWeight + treeChoice;
    dppChoice = s.dppWeight + stationaryChoice;
}

void Mcmc::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    model->regenerateLikelihood();
    model->accept();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
            std::cout << "Accept Rates Since Last Tuning Iteration:\t" << 
                         "Tree Rate=" << (double)tree->treeAcceptCount/(double)tree->treeCount << 
                         "\tBranch Rate=" << (double)tree->branchAcceptCount/(double)tree->branchCount <<
                         "\tStationary Rate=" << (double)codonMatrix->stationaryAcceptCount/(double)codonMatrix->stationaryCount <<
                         "\tK Rate=" << (double)codonMatrix->kAcceptCount/(double)codonMatrix->kCount <<
                         "\tR Rate=" << (double)codonMatrix->rAcceptCount/(double)codonMatrix->rCount <<
                         "\tOmega Rate=" << (double)dpp->omegaAcceptCount/(double)dpp->omegaCount << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        std::function<double()> updater;


        double randomMove = rng.uniformRv() * dppChoice;
        if(randomMove < kChoice){
            updater = [this]() { return codonMatrix->updateK(); };
        }
        else if(randomMove < rChoice){
            updater = [this]() { return codonMatrix->updateR(); };
        }
        else if(randomMove < omegaChoice){
            updater = [this]() { return dpp->updateOmega(); };
        }
        else if(randomMove < treeChoice){
            updater = [this]() { return tree->update(); };
        }
        else if(randomMove < stationaryChoice){
            updater = [this]() { return codonMatrix->updateStationary(); };
        }
        else if(randomMove < dppChoice){
            updater = [this]() { return dpp->updateDPP(); };
        }
        

        for(int i = 0; i < metropolisUpdates; i++){
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

        std::function<double()> updater;

        double randomMove = rng.uniformRv() * dppChoice;
        if(randomMove < kChoice){
            updater = [this]() { return codonMatrix->updateK(); };
        }
        else if(randomMove < rChoice){
            updater = [this]() { return codonMatrix->updateR(); };
        }
        else if(randomMove < omegaChoice){
            updater = [this]() { return dpp->updateOmega(); };
        }
        else if(randomMove < treeChoice){
            updater = [this]() { return tree->update(); };
        }
        else if(randomMove < stationaryChoice){
            updater = [this]() { return codonMatrix->updateStationary(); };
        }
        else if(randomMove < dppChoice){
            updater = [this]() { return dpp->updateDPP(); };
        }
        

        for(int i = 0; i < metropolisUpdates; i++){
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
    }
}
