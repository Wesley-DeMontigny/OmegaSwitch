#include "MCMC.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/Settings.hpp"
#include "Move.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

#ifdef TIME_PROFILE
#include <chrono>
#endif

/**
 * @brief Construct a new MCMC::MCMC object
 * 
 * @param m The model to sample from
 * @param mv The collection of Metropolis-Hastings proposals to make on the parameters
 * @param s The user settings
 * @param dBO Whether Bayesian optimization is being used or not
 */
MCMC::MCMC(Model* m,  std::vector<Move>& mv, Settings& s, bool dBO) : 
    model(m), optim(6, s.bayesOptFrequency), disableBayesOpt(dBO), moves(mv), totalWeight(0.0) { 
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

    for(Move& move : moves){
        totalWeight += move.weight;
    }

    model->regenerateLikelihood();
    model->accept();
}

/**
 * @brief Perform a single Gibbs iteration using a randomly chosen move
 */
double MCMC::GibbsIteration(double currentLnPosterior){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double adjustedTotal = totalWeight;
    for(Move& mv : moves){
        if(mv.condition && !mv.condition()){ // If our moves are conditional we don't want to include them in the total weight
            adjustedTotal -= mv.weight;
        }
    }

    double randomMove = rng.uniformRv() * adjustedTotal;
    std::function<double()> updater;
    double cumSum = 0.0;
    int numUpdates = 1;
    for(Move& mv : moves){
        if(mv.condition()){
            cumSum += mv.weight;
            if(randomMove < cumSum){
                updater = mv.action;
                numUpdates = mv.gibbsIterations;
                break;
            }
        }
    }


    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point initTime = std::chrono::steady_clock::now();
    #endif


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

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point finalTime = std::chrono::steady_clock::now();
    std::cout << "Gibbs iteration was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(finalTime - initTime).count() << "[milliseconds]" << std::endl;
    #endif

    return currentLnPosterior;
}

/**
 * @brief Start a burn-in run to get the Markov chain to a reasonable region 
 * of parameter space.
 */
void MCMC::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    std::vector<double> initTuning = model->getTunableParameters();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
            model->printAcceptanceRates();
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }

    if(!disableBayesOpt && bayesOptIter > 0){
        std::vector<double> currVec = model->getTunableParameters();
        std::vector<double> diffVec = initTuning;
        for(int i = 0; i < diffVec.size(); i++){
            diffVec[i] -= currVec[i];
        }
        optim.setBounds(diffVec, currVec);

        int sampleCount = 0;
        std::vector<std::vector<double>> posteriorSamples;
        for(int n = 1; n<= bayesOptIter; n++){
            currentLnPosterior = GibbsIteration(currentLnPosterior);
            std::vector<double> record = model->getTunableParameterRecord();
            posteriorSamples.push_back(record);

            if(n % bayesOptFreq == 0){
                std::cout << "Performing Bayesian Optimization..." << std::endl;
                currVec = model->getTunableParameters();
                double objective = optim.objective(posteriorSamples);
                std::cout << "Parameters had score " << objective << std::endl;
                optim.registerSample(currVec, objective);
                if(sampleCount <= 2){ // For samples 2 and 3 we just shuffle the values a bit
                    for(int i = 0; i < currVec.size(); i++){
                        currVec[i] *= std::exp((rng.uniformRv() - 0.5));
                    }
                    model->setTunableParameters(currVec);
                }
                else {
                    optim.updateGaussianProcess();
                    std::vector<double> newParams = optim.maximizeAcquisition();
                    model->setTunableParameters(newParams);
                }
                sampleCount++;
                posteriorSamples.clear();
            }
        }

        std::cout << "Choosing optimal parameters..." << std::endl;
        std::vector<double> optimParams = optim.getMaximum();
        model->setTunableParameters(optimParams);
    }
}

/**
 * @brief Start the sampling run of the Markov chain
 * 
 */
void MCMC::run(){
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    model->writeLogHeaders();

    for(int n = 1; n <= numIter; n++){
        if(n % printFreq == 0){
            model->printTabular(n);
        }
        if(n % sampleFreq == 0){
            model->writeLogData(n);
        }

        currentLnPosterior = GibbsIteration(currentLnPosterior);
    }
}
