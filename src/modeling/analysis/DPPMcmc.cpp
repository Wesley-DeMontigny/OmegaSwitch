#include "DPPMcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/DPPModel.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/DPPMatrix.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>

DPPMcmc::DPPMcmc(DPPModel* m, TreeParameter* t, DPPMatrix* cm, DirichletProcessPrior* d, Settings& s, bool dBO) : 
    model(m), dpp(d), codonMatrix(cm), tree(t), generalUpdates(3), stationaryUpdates(5), treeUpdates(0), optim(5, s.bayesOptFrequency), disableBayesOpt(dBO)  { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    
    if(!disableBayesOpt){
        bayesOptFreq = s.bayesOptFrequency;
        bayesOptIter = s.bayesOpt;
    }

    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
    dppLog = s.dppOutput;
    tipsLog = s.tipsOutput;
    ancestralLog = s.ancestralStatesOutput;
    branchLog = s.branchOutput;

    kChoice = s.kWeight;
    rChoice = s.rWeight + kChoice;
    treeChoice = s.treeWeight + rChoice;
    omegaChoice = s.omegaWeight + treeChoice;
    stationaryChoice = s.stationaryWeight + omegaChoice;
    dppChoice = s.dppWeight + stationaryChoice;

    treeUpdates = (int)(tree->getTree()->getBranchLengths().size() * 0.5);

    model->regenerateLikelihood();
    model->accept();
}

double DPPMcmc::GibbsIteration(double currentLnPosterior){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double randomMove = rng.uniformRv() * dppChoice;

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point initTime = std::chrono::steady_clock::now();
    #endif
    

    if(randomMove < stationaryChoice){
        std::function<double()> updater;
        int numUpdates = generalUpdates;

        if(randomMove < kChoice){
            updater = [this]() { return codonMatrix->updateK(); };
            #if LOGGING==1
            std::cout << "Updating K..." << std::endl;
            #endif
        }
        else if(randomMove < rChoice){
            updater = [this]() { return codonMatrix->updateR(); };
            #if LOGGING==1
            std::cout << "Updating R..." << std::endl;
            #endif
        }
        else if(randomMove < treeChoice){
            updater = [this]() { return tree->update(); };
            numUpdates = treeUpdates;
            #if LOGGING==1
            std::cout << "Updating Tree..." << std::endl;
            #endif
        }
        else if(randomMove < omegaChoice){
            updater = [this]() { return dpp->updateOmega(); };
            numUpdates = std::max(dpp->getNumCategories(), generalUpdates) * 2;
            #if LOGGING==1
            std::cout << "Updating Omega..." << std::endl;
            #endif
        }
        else if(randomMove < stationaryChoice){
            updater = [this]() { return codonMatrix->updateStationary(); };
            numUpdates = stationaryUpdates;
            #if LOGGING==1
            std::cout << "Updating Stationary..." << std::endl;
            #endif
        }

        for(int i = 0; i < numUpdates; i++){
            double lnProposalRatio = updater();
            model->regenerateLikelihood();

            double newLnPosterior = model->lnLikelihood() + model->lnPrior();

            double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
            double lnR = lnProposalRatio + lnPosteriorRatio;

            #if LOGGING==1
            std::cout << "Evalulating proposal with acceptance ratio of " << lnR << std::endl;
            #endif

            if(std::log(rng.uniformRv()) < lnR){
                #if LOGGING==1
                std::cout << "Accepted proposal!" << std::endl;
                #endif
                model->accept();
                currentLnPosterior = newLnPosterior;
            }
            else{
                #if LOGGING==1
                std::cout << "Rejected proposal!" << std::endl;
                #endif
                model->reject();
            }
        }
    }
    else {
        #if LOGGING==1
        std::cout << "Updating the DPP..." << std::endl;
        #endif
        dpp->updateDPP();

        model->regenerateLikelihood();
        currentLnPosterior = model->lnLikelihood() + model->lnPrior();
        model->accept();
    }

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point finalTime = std::chrono::steady_clock::now();
    std::cout << "Gibbs iteration was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(finalTime - initTime).count() << "[milliseconds]" << std::endl;
    #endif

    return currentLnPosterior;
}

void DPPMcmc::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    std::vector<double> initTuning(5, 0.0);
    initTuning[0] = tree->treeAlpha;
    initTuning[1] = tree->branchDelta;
    initTuning[2] = codonMatrix->stationaryAlpha;
    initTuning[3] = codonMatrix->kDelta;
    initTuning[4] = codonMatrix->rDelta;
    // We are excluding omega. That will be tuned with acceptance rate because it is really weird to handle autocorrelation in a Dirichlet Process distributed parameter


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

    if(!disableBayesOpt && bayesOptIter > 0){
        std::vector<double> diffVec = {initTuning[0] - tree->treeAlpha, initTuning[1] - tree->branchDelta, initTuning[2] - codonMatrix->stationaryAlpha, initTuning[3] - codonMatrix->kDelta, initTuning[4] - codonMatrix->rDelta};
        std::vector<double> currVec = {tree->treeAlpha, tree->branchDelta, codonMatrix->stationaryAlpha, codonMatrix->kDelta, codonMatrix->rDelta};
        optim.setBounds(diffVec, currVec);

        #if LOGGING==1
        std::cout << "Tuning parameters are:" << std::endl;
        std::cout << "\t1: " << tree->treeAlpha << std::endl;
        std::cout << "\t2: " << tree->branchDelta << std::endl;
        std::cout << "\t3: " << codonMatrix->stationaryAlpha << std::endl;
        std::cout << "\t4: " << codonMatrix->kDelta << std::endl;
        std::cout << "\t5: " << codonMatrix->rDelta << std::endl;
        #endif

        int sampleCount = 0;
        std::vector<std::vector<double>> posteriorSamples;
        for(int n = 1; n<= bayesOptIter; n++){
            currentLnPosterior = GibbsIteration(currentLnPosterior);
            std::vector<double> record = {codonMatrix->getK(), codonMatrix->getR()};
            for(double entry : codonMatrix->getRawStationary())
                record.push_back(entry);
            for(double entry : tree->getTree()->getBranchLengths())
                record.push_back(entry);

            posteriorSamples.push_back(record);

            if(n % bayesOptFreq == 0){
                std::cout << "Performing Bayesian Optimization..." << std::endl;
                currVec = {tree->treeAlpha, tree->branchDelta, codonMatrix->stationaryAlpha, codonMatrix->kDelta, codonMatrix->rDelta};
                double objective = optim.objective(posteriorSamples);
                std::cout << "Parameters had score " << objective << std::endl;
                optim.registerSample(currVec, objective);
                if(sampleCount <= 2){ // For samples 2 and 3 we just shuffle the values a bit
                    tree->treeAlpha *= std::exp(0.5* (rng.uniformRv() - 0.5));
                    tree->branchDelta *= std::exp(0.5* (rng.uniformRv() - 0.5));
                    codonMatrix->stationaryAlpha *= std::exp(0.5* (rng.uniformRv() - 0.5));
                    codonMatrix->kDelta *= std::exp(0.5* (rng.uniformRv() - 0.5));
                    codonMatrix->rDelta *= std::exp(0.5* (rng.uniformRv() - 0.5));
                }
                else {
                    optim.updateGaussianProcess();
                    std::vector<double> newParams = optim.maximizeAcquisition();
                    tree->treeAlpha = newParams[0];
                    tree->branchDelta = newParams[1];
                    codonMatrix->stationaryAlpha = newParams[2];
                    codonMatrix->kDelta = newParams[3];
                    codonMatrix->rDelta = newParams[4];
                }
                sampleCount++;
                posteriorSamples.clear();
                #if LOGGING==1
                std::cout << "Continuing MCMC..." << std::endl;
                #endif
            }
        }

        std::cout << "Choosing optimal parameters..." << std::endl;
        std::vector<double> optimParams = optim.getMaximum();
        tree->treeAlpha = optimParams[0];
        tree->branchDelta = optimParams[1];
        codonMatrix->stationaryAlpha = optimParams[2];
        codonMatrix->kDelta = optimParams[3];
        codonMatrix->rDelta = optimParams[4];
    }
}

void DPPMcmc::run(){
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

    if(branchLog != ""){
        fs.open(branchLog, std::ofstream::out);
        fs << model->branchHeader();
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

            if(branchLog != ""){
                fs.open(branchLog, std::ofstream::app);
                fs << model->branchOut(n);
                fs.close();
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
