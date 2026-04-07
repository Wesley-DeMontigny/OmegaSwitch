#include "MCMC.hpp"
#include "misc/Msg.hpp"
#include "misc/RandomVariable.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "misc/Settings.hpp"
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
MCMC::MCMC(Model* m, std::vector<Move>& mv, Settings& s, bool dBO) :
    optim(6, s.bayesOptFrequency), disableBayesOpt(dBO), hasTemperedChain(false),
    hotBeta(1.0), sampleChainIndex(0), swapAcceptCount(0), swapCount(0) {
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    swapFreq = s.mcmcmcSwapFrequency;
    models[0] = m;
    models[1] = nullptr;
    moves[0] = &mv;
    moves[1] = nullptr;
    totalWeights[0] = 0.0;
    totalWeights[1] = 0.0;
    currentLnPosteriors[0] = 0.0;
    currentLnPosteriors[1] = 0.0;

    if(!disableBayesOpt){
        bayesOptFreq = s.bayesOptFrequency;
        bayesOptIter = s.bayesOpt;
    }
    else {
        bayesOptFreq = 0;
        bayesOptIter = 0;
        numBurnIn += s.bayesOpt;
    }

    for(Move& move : *moves[0]){
        totalWeights[0] += move.weight;
    }
}

/**
 * @brief Construct a new two-chain object for MCMCMC.
 */
MCMC::MCMC(Model* coldModel, std::vector<Move>& coldMoves,
           Model* temperedModel, std::vector<Move>& temperedMoves,
           Settings& s, bool dBO) :
    optim(6, s.bayesOptFrequency), disableBayesOpt(dBO), hasTemperedChain(true),
    hotBeta(s.mcmcmcBeta), sampleChainIndex(0), swapAcceptCount(0), swapCount(0) {
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    swapFreq = s.mcmcmcSwapFrequency;
    models[0] = coldModel;
    models[1] = temperedModel;
    moves[0] = &coldMoves;
    moves[1] = &temperedMoves;
    totalWeights[0] = 0.0;
    totalWeights[1] = 0.0;
    currentLnPosteriors[0] = 0.0;
    currentLnPosteriors[1] = 0.0;
    bayesOptFreq = 0;
    bayesOptIter = 0;

    if(!disableBayesOpt && s.bayesOpt > 0){
        Msg::error("Bayesian optimization is not supported while MCMCMC is enabled.");
    }

    for(Move& move : *moves[0]){
        totalWeights[0] += move.weight;
    }
    for(Move& move : *moves[1]){
        totalWeights[1] += move.weight;
    }
}

double MCMC::initializeChain(int chainIndex){
    models[chainIndex]->regenerateLikelihood();
    models[chainIndex]->accept();
    return models[chainIndex]->lnLikelihood() + models[chainIndex]->lnPrior();
}

double MCMC::getBetaForChain(int chainIndex) const {
    if(!hasTemperedChain){
        return 1.0;
    }
    return (chainIndex == sampleChainIndex ? 1.0 : hotBeta);
}

Model* MCMC::getSampleModel() const {
    return models[sampleChainIndex];
}

Model* MCMC::getTemperedModel() const {
    return models[1 - sampleChainIndex];
}

/**
 * @brief Perform a single Gibbs iteration using a randomly chosen move
 */
double MCMC::GibbsIteration(Model* activeModel, std::vector<Move>& activeMoves,
                            double currentLnPosterior, double beta,
                            double currentTotalWeight){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    activeModel->setCountTuningEvents(beta == 1.0);

    double adjustedTotal = currentTotalWeight;
    for(Move& mv : activeMoves){
        if(!mv.condition()){
            adjustedTotal -= mv.weight;
        }
    }

    double randomMove = rng.uniformRv() * adjustedTotal;
    std::function<double()> updater;
    double cumSum = 0.0;
    int numUpdates = 1;
    for(Move& mv : activeMoves){
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
        activeModel->setPower(beta);
        double lnProposalRatio = updater();
        activeModel->regenerateLikelihood();

        double newLnPosterior = activeModel->lnLikelihood() + activeModel->lnPrior();
        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + (beta * lnPosteriorRatio);

        if(std::log(rng.uniformRv()) < lnR){
            activeModel->accept();
            currentLnPosterior = newLnPosterior;
        }
        else{
            activeModel->reject();
        }
    }

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point finalTime = std::chrono::steady_clock::now();
    std::cout << "Gibbs iteration was completed in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(finalTime - initTime).count()
              << "[milliseconds]" << std::endl;
    #endif

    return currentLnPosterior;
}

bool MCMC::attemptSwap(){
    if(!hasTemperedChain){
        return false;
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int coldIndex = sampleChainIndex;
    int hotIndex = 1 - sampleChainIndex;
    double lnR = (1.0 - hotBeta) * (currentLnPosteriors[hotIndex] - currentLnPosteriors[coldIndex]);

    swapCount += 1;
    if(std::log(rng.uniformRv()) < lnR){
        sampleChainIndex = hotIndex;
        swapAcceptCount += 1;
        return true;
    }

    return false;
}

void MCMC::printBurninAcceptanceRates() const {
    if(!hasTemperedChain){
        models[0]->printAcceptanceRates();
        return;
    }

    // THESE ARE NOT SHARED ACCEPTANCE RATES
    // What is output here are the real rates, not the rates we tune on.
    std::cout << "Posterior Move Rates: ";
    getSampleModel()->printAcceptanceRates();
    std::cout << "Tempered Move Rates: ";
    getTemperedModel()->printAcceptanceRates();
    std::cout << "Swap Acceptance Rate: "
              << (swapCount > 0 ? (double)swapAcceptCount / (double)swapCount : 0.0) << std::endl;
}

/**
 * @brief Start a burn-in run to get the Markov chain to a reasonable region
 * of parameter space.
 */
void MCMC::burnin(){
    currentLnPosteriors[0] = initializeChain(0);
    if(hasTemperedChain){
        currentLnPosteriors[1] = initializeChain(1);
    }

    std::vector<double> initTuning = getSampleModel()->getTunableParameters();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosteriors[sampleChainIndex];
            if(hasTemperedChain){
                std::cout << " | " << currentLnPosteriors[1 - sampleChainIndex];
            }
            std::cout << std::endl;
            printBurninAcceptanceRates();
        }
        if(n % tuneFreq == 0){
            models[0]->tuneMoves();
            if(hasTemperedChain){
                models[1]->tuneMoves();
            }
        }

        currentLnPosteriors[0] = GibbsIteration(models[0], *moves[0], currentLnPosteriors[0], getBetaForChain(0), totalWeights[0]);
        if(hasTemperedChain){
            currentLnPosteriors[1] = GibbsIteration(models[1], *moves[1], currentLnPosteriors[1], getBetaForChain(1), totalWeights[1]);
            if(n % swapFreq == 0){
                attemptSwap();
            }
        }
    }

    if(!disableBayesOpt && bayesOptIter > 0){
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        std::vector<double> currVec = models[0]->getTunableParameters();
        std::vector<double> diffVec = initTuning;
        for(int i = 0; i < diffVec.size(); i++){
            diffVec[i] -= currVec[i];
        }
        optim.setBounds(diffVec, currVec);

        int sampleCount = 0;
        std::vector<std::vector<double>> posteriorSamples;
        for(int n = 1; n <= bayesOptIter; n++){
            currentLnPosteriors[0] = GibbsIteration(models[0], *moves[0], currentLnPosteriors[0], 1.0, totalWeights[0]);
            std::vector<double> record = models[0]->getTunableParameterRecord();
            posteriorSamples.push_back(record);

            if(n % bayesOptFreq == 0){
                std::cout << "Performing Bayesian Optimization..." << std::endl;
                currVec = models[0]->getTunableParameters();
                double objective = optim.objective(posteriorSamples);
                std::cout << "Parameters had score " << objective << std::endl;
                optim.registerSample(currVec, objective);
                if(sampleCount <= 2){
                    for(int i = 0; i < currVec.size(); i++){
                        currVec[i] *= std::exp((rng.uniformRv() - 0.5));
                    }
                    models[0]->setTunableParameters(currVec);
                }
                else {
                    optim.updateGaussianProcess();
                    std::vector<double> newParams = optim.maximizeAcquisition();
                    models[0]->setTunableParameters(newParams);
                }
                sampleCount++;
                posteriorSamples.clear();
            }
        }

        std::cout << "Choosing optimal parameters..." << std::endl;
        std::vector<double> optimParams = optim.getMaximum();
        models[0]->setTunableParameters(optimParams);
    }
}

/**
 * @brief Start the sampling run of the Markov chain
 */
void MCMC::run(){
    currentLnPosteriors[0] = initializeChain(0);
    if(hasTemperedChain){
        currentLnPosteriors[1] = initializeChain(1);
    }

    getSampleModel()->writeLogHeaders();

    for(int n = 1; n <= numIter; n++){
        if(n % printFreq == 0){
            getSampleModel()->printTabular(n);
        }
        if(n % sampleFreq == 0){
            getSampleModel()->writeLogData(n);
        }

        currentLnPosteriors[0] = GibbsIteration(models[0], *moves[0], currentLnPosteriors[0], getBetaForChain(0), totalWeights[0]);
        if(hasTemperedChain){
            currentLnPosteriors[1] = GibbsIteration(models[1], *moves[1], currentLnPosteriors[1], getBetaForChain(1), totalWeights[1]);
            if(n % swapFreq == 0){
                attemptSwap();
            }
        }
    }
}
