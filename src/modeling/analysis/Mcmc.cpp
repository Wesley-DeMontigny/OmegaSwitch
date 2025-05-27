#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/CodonMultiMatrix.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

Mcmc::Mcmc(Model* m, TreeParameter* t, CodonMultiMatrix* cmm, DirichletProcessPrior* d, Settings& s) : 
    model(m), dpp(d), codonMatrix(cmm), tree(t), generalUpdates(5), stationaryUpdates(10), treeUpdates(0) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
    dppLog = s.dppOutput;
    tipsLog = s.tipsOutput;
    ancestralLog = s.ancestralStatesOutput;

    kChoice = s.kWeight;
    rChoice = s.rWeight + kChoice;
    treeChoice = s.treeWeight + rChoice;
    omegaChoice = s.omegaWeight + treeChoice;
    stationaryChoice = s.stationaryWeight + omegaChoice;
    dppChoice = s.dppWeight + stationaryChoice;

    treeUpdates = (int)(tree->getTree()->getBranchLengths().size() * 0.75);

    model->regenerateLikelihood();
    model->accept();
}

double Mcmc::GibbsIteration(double currentLnPosterior){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double randomMove = rng.uniformRv() * dppChoice;

    if(randomMove < stationaryChoice){
        std::function<double()> updater;
        int numUpdates = generalUpdates;

        if(randomMove < kChoice){
            updater = [this]() { return codonMatrix->updateK(); };
        }
        else if(randomMove < rChoice){
            updater = [this]() { return codonMatrix->updateR(); };
        }
        else if(randomMove < treeChoice){
            updater = [this]() { return tree->update(); };
            numUpdates = treeUpdates;
        }
        else if(randomMove < omegaChoice){
            updater = [this]() { return dpp->updateOmega(); };
            numUpdates = dpp->getNumCategories() * 2;
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
    }
    else {
        dpp->updateDPP();

        model->regenerateLikelihood();
        currentLnPosterior = model->lnLikelihood() + model->lnPrior();
        model->accept();
    }

    return currentLnPosterior;
}

void Mcmc::burnin(){
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
            std::cout << "Accept Rates Since Last Tuning Iteration:" << 
                         "\tTree Rate=" << (double)tree->treeAcceptCount/(double)tree->treeCount << 
                         "\tBranch Rate=" << (double)tree->branchAcceptCount/(double)tree->branchCount <<
                         "\tStationary Rate=" << (double)codonMatrix->stationaryAcceptCount/(double)codonMatrix->stationaryCount <<
                         "\tK Rate=" << (double)codonMatrix->kAcceptCount/(double)codonMatrix->kCount <<
                         "\tR Rate=" << (double)codonMatrix->rAcceptCount/(double)codonMatrix->rCount <<
                         "\tOmega Rate=" << (double)dpp->omegaAcceptCount/(double)dpp->omegaCount << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }
}

void Mcmc::run(){
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

    if(dppLog != ""){
        fs.open(dppLog, std::ofstream::out);
        fs << model->dppHeader();
        fs.close();
    }

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

            if(dppLog != ""){
                fs.open(dppLog, std::ofstream::app);
                fs << model->dppOut(n);
                fs.close();
                fs.clear();
            }

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
