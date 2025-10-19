/*
    EXPERIMENTAL FEATURE!
*/
#include "BayesianOptimizer.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>
#include "Math.hpp"
#include "Matrix.hpp"
#include "RandomVariable.hpp"
#include "Probability.hpp"
#include "Msg.hpp"

/**
 * @brief For each element x in v, return a matrix with elements 10^x. Used to convert between log spaces.
 */
std::vector<double> pow10Vector(std::vector<double> v){
    std::vector<double> returnVec;
    for(double n : v)
        returnVec.push_back(std::pow(10.0,n));
    return returnVec;
}

/**
 * @brief Constructor for the Bayesian Optimizer. Takes in the number of parameters we are trying to optimize and the
 * number of itersations per sample. 
 */
BayesianOptimizer::BayesianOptimizer(int nP, int s) : numParams(nP), iterationsPerSample(s), hyperparams(nP, 1.0) {}

/**
 * @brief Compute the autocorrelation score within a particular window.
 */
double BayesianOptimizer::autocorrelationScore(const std::vector<std::vector<double>>& r, int start, int end){
    int window = 25;

    std::vector<double> mean(r[0].size(), 0.0);
    for(int i = start; i < end; i++)
        for(int j = 0; j < r[0].size(); j++)
            mean[j] += r[i][j];
    for(int i = 0; i < mean.size(); i++)
        mean[i] /= window;

    std::vector<double> variance(r[0].size(), 0.0);
    for(int i = start; i < end; i++)
        for(int j = 0; j < r[0].size(); j++)
            variance[j] += std::pow(r[i][j] - mean[j], 2);
    
    double var_sum = 0.0;
    for(int i = 0; i < mean.size(); i++){
        variance[i] /= window;
        var_sum += variance[i];
    }

    double autoCorrSum = 0.0;

    for(int i = 1; i < window; i++){
        double sum = 0;
        for(int t = start; t < end - i; t++){
            for(int j = 0; j < r[0].size(); j++)
                sum += (r[t][j] - mean[j]) * (r[t+i][j] - mean[j]);
        }
        sum /= (window - i)*var_sum;

        autoCorrSum += std::abs(sum);
    }

    return 1 - (autoCorrSum / window - 1);
}

/**
 * @brief Get the set of MCMC tunable parameters with the maximum score (miniumum autocorrelation)
 */
std::vector<double> BayesianOptimizer::getMaximum(){

    std::sort(samples.begin(), samples.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
        return a.score > b.score;
    });


    std::cout << "Using parameters with maximum score of " << samples[0].score << std::endl;;

    return samples[0].params;
}

/**
 * @brief Get the top N set of MCMC tunable parameters by autocorrelation score
 */
std::vector<std::vector<double>> BayesianOptimizer::getMaximumN(int N){
    std::vector<std::vector<double>> returnVec;

    std::sort(samples.begin(), samples.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
        return a.score > b.score;
    });
    
    for(int i = 0; i < N; i++){
        returnVec.push_back(samples[i].params);
    }

    return returnVec;
}

/**
 * @brief Compute the entry in the covariance matrix corresponding to the Kernel value for
 * two tunable parameter sets.
 */
double BayesianOptimizer::kernel(std::vector<double>& a, std::vector<double>& b){
    double sum = 0.0;
    for(int i = 0; i < numParams; i++){
        double diff = a[i] - b[i];
        sum += std::pow(diff, 2) / (std::max(1e-10, std::pow(hyperparams[i], 2))); // To make sure we don't underflow
    }
    sum /= -2.0;

    return sampleVariance * std::exp(sum);
}

/**
 * @brief 
 */
double BayesianOptimizer::marginalLogLikelihood(){
    int numSamples = samples.size();
    double likelihood = -0.5;

    double sum = 0.0;
    double logDet = 0.0;
    for(int i = 0; i < numSamples; i++){
        sum += samples[i].score * alpha[i];
        logDet += std::log(choleskyFactor(i,i));
    }
    likelihood *= sum;
    likelihood -= logDet;

    likelihood -= numSamples * 0.5 * std::log(M_PI * 2);

    return likelihood;
}

/**
 * @brief Compute the gradient of the marginal log likelihood of the Gaussian process
 * with respect to the hyperparameters. See Rasmussen and Williams Ch. 5 for theoretical background.
 */
std::vector<double> BayesianOptimizer::nLLGradient(){
    int numSamples = samples.size();
    std::vector<double> gradients(numParams, 0.0);

    std::vector<std::vector<double>> kernelInvColumns;

    Matrix<double> choleskyFactorTranspose(numSamples, numSamples, 0.0);
    Math::transposeMatrix(choleskyFactor, choleskyFactorTranspose);

    for(int i = 0; i < numSamples; i++){
        kernelInvColumns.push_back(std::vector<double>(numSamples, 0.0));
        kernelInvColumns[i][i] = 1.0;
        Math::forwardSubstitutionRow(choleskyFactor, kernelInvColumns[i]);

        Math::backSubstitutionRow(choleskyFactorTranspose, kernelInvColumns[i]);
    }

    Matrix<double> difference(numSamples, numSamples, 0.0);
    for(int i = 0; i < numSamples; i++){
        for(int j = 0; j < numSamples; j++){
            difference(i,j) = alpha[i] * alpha[j];
        }
    }

    for(int i = 0; i < numSamples; i++){
        for(int j = 0; j < numSamples; j++){
            difference(i,j) -= kernelInvColumns[j][i];
        }
    }

    // Gradient of the Kernel with respect to the hyperparameters
    for(int p = 0; p < numParams; p++){
        Matrix<double> kCopy = kernelMatrix.copy();

        double hyperParamCubed = std::pow(hyperparams[p], 3.0);

        for(int i = 0; i < numSamples; i++){
            for(int j = 0; j < numSamples; j++){
                double diff = samples[i].params[p] - samples[j].params[p];
                kCopy(i,j) *= std::pow(diff, 2.0);
            }
        }

        kCopy /= hyperParamCubed;

        Matrix<double> matProd = difference * kCopy;
        double trace = 0.0;
        for(int i = 0; i < numSamples; i++){
            trace += matProd(i,i);
        }

        gradients[p] = trace * -0.5; // Multiply by -1.0 this for the negative log likelihood
    }
    
    return gradients;
}

/**
 * @brief Compute the log prior probability of the hyperparameters. Here we are
 * just using a standard log Normal
 */
double BayesianOptimizer::logPrior(){
    double total = 0.0;
    for(int i = 0; i < numParams; i++){
        double priorProbs = log(hyperparams[i]) + 
                            std::log(std::sqrt(2.0 * M_PI)) +
                            std::pow(std::log(hyperparams[i]) - 0, 2.0)/2.0;
        total -= priorProbs;
    }

    return total;
}

/**
 * @brief Compute the gradient of the log prior with respect to the hyperparameter
 * values.
 */
std::vector<double> BayesianOptimizer::nLPGradient(){
    std::vector<double> gradients(numParams, 0.0);

    // Gradients for the log-normal prior
    for(int i = 0; i < numParams; i++){
        double grad = 1/hyperparams[i];
        grad *= 1.0 + ((std::log(hyperparams[i]) - 0.0)); // We are putting a regularizing prior on this so mean of 0 and variance of 1
        gradients[i] = grad;
    }

    return gradients;
}

/**
 * @brief Register a set of tunable parameters (s) and their corresponding
 * value in the GP (the autocorrelation score).
 */
void BayesianOptimizer::registerSample(std::vector<double> s, double o){
    ParamScorePair sample;
    sample.score = o;
    sample.params = s;
    samples.push_back(sample);

    sampleMean = 0.0;
    for(int i = 0; i < samples.size(); i++){
        sampleMean += samples[i].score;
    }
    sampleMean /= samples.size();

    sampleVariance = 0.0;
    for(int i = 0; i < samples.size(); i++){
        sampleVariance += pow(samples[i].score - sampleMean, 2);
    }
    sampleVariance /= samples.size();

    meanSampleComponents.clear();
    for(int i = 0; i < numParams; i++){
        meanSampleComponents.push_back(0.0);

        for(int j = 0; j < samples.size(); j++)
            meanSampleComponents[i] += samples[j].params[i];
        
        meanSampleComponents[i] /= samples.size();
    }
}

/**
 * @brief Compute the upper confidence bound of a set of hypothetical tunable parameters
 * using an exploration parameter of 0.25.
 */
double BayesianOptimizer::UCB(std::vector<double>& sample, int numSamples){
    std::vector<double> kStar(numSamples, 0.0);
    for(int j = 0; j < numSamples; j++){
        kStar[j] = kernel(sample, samples[j].params);
    }

    double predictive_mean = 0.0;
    for(int j = 0; j < numSamples; j++){
        predictive_mean += kStar[j] * alpha[j];
    }

    std::vector<double> v = kStar;
    Math::forwardSubstitutionRow(choleskyFactor, v);

    double predictive_variance = kernel(kStar, kStar);
    for(int j = 0; j < numSamples; j++){
        predictive_variance -= v[j]*v[j];
    }
    double predictive_std = std::sqrt(predictive_variance);

    // Set exploration parameter to 0.25
    return predictive_mean + 0.25*predictive_std;
}

/**
 * @brief Select the best next set of tunable MCMC parameters based on the Gaussian Process
 * and the UCB acquisition function. Right now we are just looking acroos the whole space 
 * using latin hypercube sampling.
 */
std::vector<double> BayesianOptimizer::maximizeAcquisition(){
    assert(samples.size() >= 3); // We need to make sure that we have selected three coherent points

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    #if LOGGING==1
    std::cout << "Using latin hypercube sampling to sample candidate parameters..." << std::endl;
    #endif

    int numLHCSamples = 2500;
    int numSamples = samples.size();

    std::vector<std::vector<double>> lhcSamples;
    std::vector<int> randomIndices;
    for(int i = 0; i < numLHCSamples; i++){
        lhcSamples.push_back(std::vector<double>(numParams, 0.0));
        randomIndices.push_back(i);
    }

    // Do Latin Hypercube Sampling
    for(int d = 0; d < numParams; d++){
        std::vector<int> availableDraws = randomIndices;
        double lB = lowerBound[d];
        double increment = (upperBound[d] - lB)/(double)(numLHCSamples);
        for(int n = 0; n < numLHCSamples; n++){
            int assigningIndex = availableDraws[(int)(rng.uniformRv() * availableDraws.size())];
            auto it = std::find(availableDraws.begin(), availableDraws.end(), assigningIndex);
            availableDraws.erase(it);

            lhcSamples[assigningIndex][d] = lB + increment*n + (rng.uniformRv() * increment);
        }
    }

    std::sort(samples.begin(), samples.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
        return a.score > b.score;
    });

    double currentMaxObjective = samples[0].score;

    #if LOGGING==1
    std::cout << "Maximizing acquisition function..." << std::endl;
    #endif

    std::vector<ParamScorePair> parameterScorePairs;
    for(int i = 0; i < numLHCSamples; i++){
        
        double score = UCB(lhcSamples[i], numSamples);

        ParamScorePair point;
        point.params = lhcSamples[i];
        point.score = score;
        parameterScorePairs.push_back(point);
    }

    std::sort(parameterScorePairs.begin(), parameterScorePairs.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
        return a.score > b.score;
    });

    return parameterScorePairs[0].params;
}

/**
 * @brief  Define the bounds of our search using the change of the parameters so far and the current center
 */
void BayesianOptimizer::setBounds(std::vector<double>& diff, std::vector<double>& mean){
    assert(diff.size() == mean.size());
    upperBound.clear();
    lowerBound.clear();

    for(int i = 0; i < diff.size(); i++){
        double abs_diff = std::min(std::abs(diff[i]), 0.75*mean[i]);
        upperBound.push_back(mean[i] + abs_diff);
        lowerBound.push_back(std::max(mean[i] - abs_diff, 1e-3));
    }

    std::cout << "Setting bounds on the MCMC parameter optimization:" << std::endl;
    for(int i = 0; i < upperBound.size(); i++){
        std::cout << "\t" << i << ": " << lowerBound[i] << "-" << upperBound[i] << std::endl;
    }
}

/**
 * @brief Update the Cholesky factor for the the covariance matrix of
 * the covariance matrix of the Gaussian process. See Rasmussen and Williams (2008) 
 * Algorithm 2.1 for details.
 */ 
 void BayesianOptimizer::updateCholesky(){
    int numSamples = samples.size();

    Matrix<double> K(numSamples, numSamples, 0.0);

    double min = INFINITY;
    for(int i = 0; i < numSamples; i++){
        for(int j = 0; j < numSamples; j++){
            double entry = kernel(samples[i].params, samples[j].params);
            K(i, j) = entry;
            if(i == j && entry < min)
                min = entry;
        }
    }

    kernelMatrix = K.copy();

    for(int i = 0; i < numSamples; i++){
        // The jitter should be on the same scale as the diagonal
        K(i, i) += 1e-6 * min;
    }

    choleskyFactor = Matrix<double>(numSamples, numSamples, 0.0);
    Math::choleskyDecomposition(K, choleskyFactor);

    std::vector<double> objectives;
    for(ParamScorePair s : samples)
        objectives.push_back(s.score);
    alpha = objectives; // This will be overwritten
    Math::forwardSubstitutionRow(choleskyFactor, alpha);

    Matrix<double> LT(numSamples, numSamples, 0.0);
    int tSuccess = Math::transposeMatrix(choleskyFactor, LT);
    Math::backSubstitutionRow(LT, alpha);
}

/**
 * @brief Update the hyperparameters of the Gaussian process using latin hypercube sampling as initialization
 * and BFGS to get to an optimum. 
 */
void BayesianOptimizer::updateGaussianProcess() {
    assert(samples.size() >= 3); // We need to make sure that we have selected three coherent points

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // Initialize simplex using latin hypercube sampling
    #if LOGGING==1
    std::cout << "Using latin hypercube sampling to get initial hyperparameter candidates..." << std::endl;
    #endif
    std::vector<ParamScorePair> hyperparamPoints;

    int numLHCSamples = 1000;
    int numSamples = samples.size();

    std::vector<std::vector<double>> lhcSamples;
    std::vector<int> randomIndices;
    for(int i = 0; i < numLHCSamples; i++){
        lhcSamples.push_back(std::vector<double>(numParams, 0.0));
        randomIndices.push_back(i);
    }

    for(int d = 0; d < numParams; d++){
        std::vector<int> availableDraws = randomIndices;
        double lB = std::floor(std::log10(lowerBound[d])) - 2.0;
        double increment = std::floor(std::log10(upperBound[d]) - lB)/(double)numLHCSamples;
        for(int n = 0; n < numLHCSamples; n++){
            int assigningIndex = availableDraws[(int)(rng.uniformRv() * availableDraws.size())];
            auto it = std::find(availableDraws.begin(), availableDraws.end(), assigningIndex);
            availableDraws.erase(it);

            lhcSamples[assigningIndex][d] = lB + increment*n + (rng.uniformRv() * increment);
        }
    }

    for(int i = 0; i < numLHCSamples; i++){
        hyperparams = pow10Vector(lhcSamples[i]);
        updateCholesky();
        double logL = -1.0*marginalLogLikelihood() - logPrior();
        ParamScorePair newVertex;
        newVertex.params = lhcSamples[i];
        newVertex.score = logL;
        hyperparamPoints.push_back(newVertex);
    }

    ParamScorePair currentValue;
    currentValue.score = INFINITY;
    for(ParamScorePair psp : hyperparamPoints){ // Minimize negative log likelihood
        if(psp < currentValue)
            currentValue = psp;
    }
    hyperparams = pow10Vector(currentValue.params);
    updateCholesky();

    currentValue.params = pow10Vector(currentValue.params);

    std::vector<std::vector<double>> history;
    std::vector<std::vector<double>> grads;
    int history_length = 5;
    double sufficientDecrease = 1e-4;
    int numIter = 0;
    bool converged = false;

    history.push_back(currentValue.params);
    grads.push_back(nLLGradient());

    std::vector<double> priorGrad = nLPGradient();
    for(int i = 0; i < numParams; i++){
        grads[0][i] += priorGrad[i];
    }

    do {
        std::vector<double> q = grads.back();;

        if(history.size() >= 2){
            std::vector<double> rhoVec;
            std::vector<double> alphaVec;

            for(int i = 1; i < grads.size(); i++){
                double rhoI = 0.0;
                for(int j = 0; j < numParams; j++){
                    rhoI += (grads[i][j] - grads[i-1][j]) * (history[i][j] - history[i-1][j]);
                }
                rhoI = 1/rhoI;
                rhoVec.push_back(rhoI);

                double alphaI = 0.0;
                for(int j = 0; j < numParams; j++){
                    alphaI += (history[i][j] - history[i-1][j]) * q[j];
                }
                alphaI *= rhoI;
                alphaVec.push_back(alphaI);

                for(int j = 0; j < numParams; j++){
                    q[j] = q[j] - alphaI * (grads[i][j] - grads[i-1][j]);
                }
            }

            double numerator = 0.0;
            double denominator = 0.0;
            for(int i = 0; i < numParams; i++){
                double y = (grads[grads.size() - 1][i] - grads[grads.size() - 2][i]);
                double s = (history[grads.size() - 1][i] - history[grads.size() - 2][i]);
                numerator += y * s;
                denominator += y * y;
            }
            double gamma = numerator/(denominator + 1e-8);
            gamma = std::clamp(gamma, 1e-4, 1e2); // For stability

            for(int i = 0; i < numParams; i++){
                q[i] *= gamma;
            }

            int k = grads.size() - 1;  // Most recent index

            for (int i = 0; i < k; i++) {
                std::vector<double> s(numParams), y(numParams);
                for (int j = 0; j < numParams; j++) {
                    s[j] = history[i+1][j] - history[i][j];
                    y[j] = grads[i+1][j] - grads[i][j];
                }

                double rho = rhoVec[i];
                double alpha = alphaVec[i];

                double beta = 0.0;
                for (int j = 0; j < numParams; j++) {
                    beta += y[j] * q[j];
                }
                beta *= rho;

                for (int j = 0; j < numParams; j++) {
                    q[j] += s[j] * (alpha - beta);
                }
            }
        }

        #if LOGGING==1
        std::cout << "Step Vector: " << std::endl;
        for(int i = 0; i < numParams; i++){
            std::cout << "\t" << i << ": " << q[i] << std::endl;
        }
        #endif

        //q now contains a valid step in the negative direction. Now we do a Line Search to find the step size that will satisfy Armijo's condition
        double epsilon = 1.0;
        double tempLogL = currentValue.score;
        std::vector<double> newHyperparams;
        bool success = true;
        while (true) {
            newHyperparams = currentValue.params;
            for(int i = 0; i < numParams; i++){
                newHyperparams[i] -= epsilon * q[i];
            }
            hyperparams = newHyperparams;
            updateCholesky();
            tempLogL = -1.0*marginalLogLikelihood() - logPrior();

            double condition = 0.0;
            for(int i = 0; i < numParams; i++){
                condition += grads.back()[i] * q[i];
            }

            condition *= 1.0 * sufficientDecrease * epsilon;
            condition += currentValue.score;

            if (tempLogL <= condition) {
                break;  // Armijo condition met
            }
            
            epsilon *= 0.5;  // Shrink step size

            if (epsilon < 1e-10) {
                #if LOGGING == 1
                std::cout << "Warning: Line search failed to find sufficient decrease!" << std::endl;
                #endif
                success = false;
                break;
            }
        }

        if(success == false){
            newHyperparams = currentValue.params;
            for(int i = 0; i < numParams; i++){
                newHyperparams[i] -= q[i] * 1e-10;
            }
            hyperparams = newHyperparams;
            updateCholesky();
            tempLogL = -1.0*marginalLogLikelihood() - logPrior();
            if(tempLogL >= currentValue.score){
                #if LOGGING == 1
                std::cout << "Rejecting L-BFGS Step!" << std::endl;
                #endif
                newHyperparams = currentValue.params;
                hyperparams = newHyperparams;
                updateCholesky();
                tempLogL = currentValue.score;
            }
        }

        #if LOGGING==1
        std::cout << "Step-taken (" << tempLogL << "): " << std::endl;
        for(int i = 0; i < numParams; i++){
            std::cout << "\t" << i << ": " << newHyperparams[i] - currentValue.params[i] << std::endl;
        }
        #endif

        currentValue.score = tempLogL;
        currentValue.params = newHyperparams;
        
        history.push_back(newHyperparams);
        grads.push_back(nLLGradient());

        std::vector<double> pG = nLPGradient();
        for(int i = 0; i < numParams; i++){
            grads[grads.size()-1][i] += pG[i];
        }

        // Only maintain a certain length of entries
        if(history.size() > history_length){
            history.erase(history.begin());
            grads.erase(grads.begin());
        }
        // If norm of difference is less than 1e-4 we quit
        if(history.size() > 1){
            double sum = 0.0;
            for(int i = 0; i < numParams; i++){
                sum += std::pow(history.back()[i] - history[history.size()-2][i], 2);
            }
            if(std::sqrt(sum) < 1e-4){
                converged = true;
            }
        }
        if(numIter > 10000){
            std::cout << "Warning: L-BFGS failed to converge!" << std::endl;
            break;
        }
        numIter++;
    }
    while(converged == false);

    std::cout << "Optimization has finished in " << numIter << " iterations with parameters:" << std::endl;
    for(int i = 0; i < numParams; i++){
        std::cout << "\t" << i << ": " << hyperparams[i] << std::endl;
    }
}

/**
 * @brief Smooth the autocorrelation score over the entire MCMC trajectory.
 */
double BayesianOptimizer::smoothAverageAutocorrelation(const std::vector<std::vector<double>>& r){
    int L = r.size();

    assert(L >= 25);

    double sum = 0;
    for(int i = 25; i < L; i++){
        sum += autocorrelationScore(r, i-25, i);
    }

    return sum / ((double)(L) - 25.0 + 1.0);
}