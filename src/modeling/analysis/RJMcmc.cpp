#include "RJMcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/RJModel.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/RJMatrix.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

RJMcmc::RJMcmc(RJModel* m, TreeParameter* t, RJMatrix* cm, Settings& s, bool dBO) : 
    model(m), codonMatrix(cm), tree(t), generalUpdates(3), stationaryUpdates(5), treeUpdates(0), optim(6, s.bayesOptFrequency), disableBayesOpt(dBO) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;

    if(!disableBayesOpt){
        bayesOptFreq = s.bayesOptFrequency;
        bayesOptIter = s.bayesOpt;
    }
    else {
        numBurnIn += s.bayesOpt;
    }

    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
    tipsLog = s.tipsOutput;
    ancestralLog = s.ancestralStatesOutput;
    branchLog = s.branchOutput;

    kChoice = s.kWeight;
    treeChoice = s.treeWeight + kChoice;
    omegaChoice = s.omegaWeight + treeChoice;
    rChoice = s.rWeight + omegaChoice;
    stationaryChoice = s.stationaryWeight + rChoice;
    rjChoice += stationaryChoice + s.rjWeight;


    treeUpdates = (int)(tree->getTree()->getBranchLengths().size() * 0.5);

    model->regenerateLikelihood();
    model->accept();
}

double RJMcmc::GibbsIteration(double currentLnPosterior){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double randomMove = rng.uniformRv() * rjChoice;

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point initTime = std::chrono::steady_clock::now();
    #endif

    std::function<double()> updater;
    int numUpdates = generalUpdates;

    if(randomMove < kChoice){
        updater = [this]() { return codonMatrix->updateK(); };
        #if LOGGING==1
        std::cout << "Updating K..." << std::endl;
        #endif
    }
    else if(randomMove < treeChoice){
        updater = [this]() { return tree->update(); };
        numUpdates = treeUpdates;
        #if LOGGING==1
        std::cout << "Updating Tree..." << std::endl;
        #endif
    }
    else if(randomMove < omegaChoice || (randomMove < rChoice && codonMatrix->getActiveOmegas() == 1)){
        updater = [this]() { return codonMatrix->updateOmega(); };
        numUpdates = codonMatrix->getActiveOmegas()*generalUpdates;
        #if LOGGING==1
        std::cout << "Updating Omega..." << std::endl;
        #endif
    }
    else if(randomMove < rChoice){
        updater = [this]() { return codonMatrix->updateR(); };
        numUpdates = generalUpdates;
        #if LOGGING==1
        std::cout << "Updating R..." << std::endl;
        #endif
    }
    else if(randomMove < stationaryChoice){
        updater = [this]() { return codonMatrix->updateStationary(); };
        numUpdates = stationaryUpdates;
        #if LOGGING==1
        std::cout << "Updating Stationary..." << std::endl;
        #endif
    }
    else if(randomMove < rjChoice){
        updater = [this]() { return codonMatrix->updateActiveOmegas(); };
        numUpdates = 1;
        #if LOGGING==1
        std::cout << "Updating Model..." << std::endl;
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

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point finalTime = std::chrono::steady_clock::now();
    std::cout << "Gibbs iteration was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(finalTime - initTime).count() << "[milliseconds]" << std::endl;
    #endif

    return currentLnPosterior;
}

void RJMcmc::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    std::vector<double> initTuning(6, 0.0);
    initTuning[0] = tree->treeAlpha;
    initTuning[1] = tree->branchDelta;
    initTuning[2] = codonMatrix->stationaryAlpha;
    initTuning[3] = codonMatrix->kDelta;
    initTuning[4] = codonMatrix->omegaDelta;
    initTuning[5] = codonMatrix->rDelta;

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
            std::cout << "Accept Rates Since Last Tuning Iteration:" << 
                         "\tTree Rate=" << (double)tree->treeAcceptCount/(double)tree->treeCount << 
                         "\tBranch Rate=" << (double)tree->branchAcceptCount/(double)tree->branchCount <<
                         "\tStationary Rate=" << (double)codonMatrix->stationaryAcceptCount/(double)codonMatrix->stationaryCount <<
                         "\tK Rate=" << (double)codonMatrix->kAcceptCount/(double)codonMatrix->kCount <<
                         "\tOmega Rate=" << (double)codonMatrix->omegaAcceptCount/(double)codonMatrix->omegaCount <<
                         "\tR Rate=" << (double)codonMatrix->rAcceptCount/(double)codonMatrix->rCount << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }

    if(!disableBayesOpt && bayesOptIter > 0){
        std::vector<double> diffVec = {initTuning[0] - tree->treeAlpha, initTuning[1] - tree->branchDelta, initTuning[2] - codonMatrix->stationaryAlpha, initTuning[3] - codonMatrix->kDelta, initTuning[4] - codonMatrix->omegaDelta, initTuning[5] - codonMatrix->rDelta};
        std::vector<double> currVec = {tree->treeAlpha, tree->branchDelta, codonMatrix->stationaryAlpha, codonMatrix->kDelta, codonMatrix->omegaDelta, codonMatrix->rDelta};
        optim.setBounds(diffVec, currVec);

        #if LOGGING==1
        std::cout << "Tuning parameters are:" << std::endl;
        std::cout << "\t1: " << tree->treeAlpha << std::endl;
        std::cout << "\t2: " << tree->branchDelta << std::endl;
        std::cout << "\t3: " << codonMatrix->stationaryAlpha << std::endl;
        std::cout << "\t4: " << codonMatrix->kDelta << std::endl;
        std::cout << "\t5: " << codonMatrix->omegaDelta << std::endl;
        std::cout << "\t6: " << codonMatrix->rDelta << std::endl;
        #endif

        int sampleCount = 0;
        std::vector<std::vector<double>> posteriorSamples;
        for(int n = 1; n<= bayesOptIter; n++){
            currentLnPosterior = GibbsIteration(currentLnPosterior);
            std::vector<double> record = {codonMatrix->getK(), codonMatrix->getOmega1(), codonMatrix->getOmega2(), codonMatrix->getOmega3(), codonMatrix->getR()};
            for(double entry : codonMatrix->getRawStationary())
                record.push_back(entry);
            for(double entry : tree->getTree()->getBranchLengths())
                record.push_back(entry);

            posteriorSamples.push_back(record);

            if(n % bayesOptFreq == 0){
                std::cout << "Performing Bayesian Optimization..." << std::endl;
                currVec = {tree->treeAlpha, tree->branchDelta, codonMatrix->stationaryAlpha, codonMatrix->kDelta, codonMatrix->omegaDelta, codonMatrix->rDelta};
                double objective = optim.objective(posteriorSamples);
                std::cout << "Parameters had score " << objective << std::endl;
                optim.registerSample(currVec, objective);
                if(sampleCount <= 2){ // For samples 2 and 3 we just shuffle the values a bit
                    tree->treeAlpha *= std::exp((rng.uniformRv() - 0.5));
                    tree->branchDelta *= std::exp((rng.uniformRv() - 0.5));
                    codonMatrix->stationaryAlpha *= std::exp((rng.uniformRv() - 0.5));
                    codonMatrix->kDelta *= std::exp((rng.uniformRv() - 0.5));
                    codonMatrix->omegaDelta *= std::exp((rng.uniformRv() - 0.5));
                    codonMatrix->rDelta *= std::exp((rng.uniformRv() - 0.5));
                }
                else {
                    optim.updateGaussianProcess();
                    std::vector<double> newParams = optim.maximizeAcquisition();
                    tree->treeAlpha = newParams[0];
                    tree->branchDelta = newParams[1];
                    codonMatrix->stationaryAlpha = newParams[2];
                    codonMatrix->kDelta = newParams[3];
                    codonMatrix->omegaDelta = newParams[4];
                    codonMatrix->rDelta = newParams[5];
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
        codonMatrix->omegaDelta = optimParams[4];
        codonMatrix->rDelta = optimParams[5];
    }
}

void RJMcmc::run(){
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
