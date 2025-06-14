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

std::vector<double> pow10Vector(std::vector<double> v){
    std::vector<double> returnVec;
    for(double n : v)
        returnVec.push_back(std::pow(10.0, n));
    return returnVec;
}

BayesianOptimizer::BayesianOptimizer(int nP, int s) : numParams(nP), iterationsPerSample(s), hyperparams(nP, 1.0) {}

BayesianOptimizer::~BayesianOptimizer() {}

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

double BayesianOptimizer::averageProportionalJumpingDistance(const std::vector<std::vector<double>>& r){
    double apjd = 0.0; // Average proportional jumping distance
    for(int i = 1; i < r.size(); i++){
        double sum = 0.0;
        for(int j = 0; j < r[0].size(); j++){
            sum += std::abs(r[i][j] - r[i-1][j])/r[i-1][j];
        }
        apjd += sum;
    }

    return apjd/r.size();
}

std::vector<double> BayesianOptimizer::getMaximum(){

    std::sort(samples.begin(), samples.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
        return a.score > b.score;
    });


    std::cout << "Using parameters with maximum score of " << samples[0].score << std::endl;;

    return samples[0].params;
}

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

double BayesianOptimizer::kernel(std::vector<double>& a, std::vector<double>& b){
    double sum = 0.0;
    for(int i = 0; i < numParams; i++){
        double diff = a[i] - b[i];
        sum += std::pow(diff, 2) / std::pow(hyperparams[i], 2);
    }
    sum /= -2.0;

    return sampleVariance * std::exp(sum);
}

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

// See Rasmussen and Williams Ch. 5
std::vector<double> BayesianOptimizer::marginalLogLikelihoodGradient(){
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

    for(int p = 0; p < numParams; p++){
        Matrix<double> kCopy = kernelMatrix.copy();
        double hyperParamCubed = std::pow(hyperparams[p], 3);

        for(int i = 0; i < numSamples; i++){
            for(int j = 0; j < numSamples; j++){
                double diff = samples[i].params[p] - samples[j].params[p];
                kCopy(i,j) *= std::pow(diff, 2);
            }
        }

        kCopy /= hyperParamCubed;

        Matrix<double> matProd = difference * kCopy;
        double trace = 0.0;
        for(int i = 0; i < numSamples; i++){
            trace += matProd(i,i);
        }

        gradients[p] = trace * 0.5;
    }
    
    return gradients;
}

void BayesianOptimizer::registerSample(std::vector<double> s, double o){
    ParamScorePair sample;
    sample.score = o;
    sample.params = s;
    samples.push_back(sample); 

    double sampleMean = 0.0;
    for(int i = 0; i < samples.size(); i++){
        sampleMean += samples[i].score;
    }
    sampleMean /= samples.size();

    sampleVariance = 0.0;
    for(int i = 0; i < samples.size(); i++){
        sampleVariance += pow(samples[i].score - sampleMean, 2);
    }
    sampleVariance /= samples.size();
}

double BayesianOptimizer::expectedImprovement(std::vector<double>& sample, double currentMaxObjective, int numSamples){
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

    // Compute expected improvement and compare
    double diff = predictive_mean - currentMaxObjective;
    double z = diff/predictive_std;
    double EI = std::max(0.0, diff) + 
                predictive_std * Probability::Normal::pdf(0, 1, z) - 
                std::abs(diff) * Probability::Normal::cdf(0, 1, z);
    
    return EI;
}

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

    // Set kappa to 0.5? We are encouraging it to stay in a somewhat reasonable area
    return predictive_mean + 0.5*predictive_std;
}

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
        
        //double score = expectedImprovement(lhcSamples[i], currentMaxObjective, numSamples);
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

void BayesianOptimizer::setBounds(std::vector<double>& diff, std::vector<double>& mean){
    assert(diff.size() == mean.size());
    upperBound.clear();
    lowerBound.clear();

    for(int i = 0; i < diff.size(); i++){
        double abs_diff = std::max(std::abs(diff[i]), 0.75*mean[i]);
        upperBound.push_back(mean[i] + abs_diff);
        lowerBound.push_back(std::max(mean[i] - abs_diff, 1e-3));
    }

    std::cout << "Setting bounds on the MCMC parameter optimization:" << std::endl;
    for(int i = 0; i < upperBound.size(); i++){
        std::cout << "\t" << i << ": " << lowerBound[i] << "-" << upperBound[i] << std::endl;
    }
}

// See Rasmussen and Williams (2008) Algorithm 2.1
void BayesianOptimizer::updateCholesky(){
    int numSamples = samples.size();

    Matrix<double> K(numSamples, numSamples, 0.0);
    for(int i = 0; i < numSamples; i++){
        for(int j = 0; j < numSamples; j++){
            K(i, j) = kernel(samples[i].params, samples[j].params);
            if(i == j)
                K(i, j) += 0.01; // For numerical stability
        }
    }

    kernelMatrix = K.copy();

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
        double lB = std::max(std::floor(std::log10(lowerBound[d]))-2.0, -2.0);
        double increment = std::max(std::floor(std::log10(upperBound[d]) - 1.0 - lB), 2.0)/(double)numLHCSamples;
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
        double marginalL = marginalLogLikelihood();
        ParamScorePair newVertex;
        newVertex.params = lhcSamples[i];
        newVertex.score = marginalL;
        hyperparamPoints.push_back(newVertex);
    }

    ParamScorePair currentValue;
    for(ParamScorePair psp : hyperparamPoints){
        if(currentValue < psp)
            currentValue = psp;
    }
    currentValue.params = pow10Vector(currentValue.params);
    hyperparams = currentValue.params;
    updateCholesky();

    std::vector<std::vector<double>> history;
    std::vector<std::vector<double>> grads;
    int history_length = 10;
    double sufficientIncrease = 1e-4;
    int numIter = 0;
    bool converged = false;

    history.push_back(currentValue.params);
    grads.push_back(marginalLogLikelihoodGradient());
    double currentLogL = marginalLogLikelihood();

    do {
        std::vector<double> q = grads.back();;

        std::vector<double> rhoVec;
        std::vector<double> alphaVec;

        for(int i = grads.size()-1; i >= 1; i--){
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

        if(history.size() > 2){
            double numerator = 0.0;
            double denominator = 0.0;
            for(int i = 0; i < numParams; i++){
                double y = (grads[grads.size() - 2][i] - grads[grads.size() - 3][i]);
                double s = (history[grads.size() - 2][i] - history[grads.size() - 3][i]);
                numerator += y * s;
                denominator += y * y;
            }
            double gamma = numerator/denominator;

            for(int i = 0; i < numParams; i++){
                q[i] *= gamma;
            }

            int k = grads.size() - 1;  // Most recent index

            for (int i = 1; i <= k; i++) {
                int offset = k - i; // The indexing is weird here because of how I push it

                double beta = 0.0;
                for (int j = 0; j < numParams; j++) {
                    beta += (grads[i][j] - grads[i-1][j]) * q[j];
                }
                beta *= rhoVec[offset];

                for (int j = 0; j < numParams; j++) {
                    q[j] += (history[i][j] - history[i-1][j]) * (alphaVec[offset] - beta);
                }
            }
        }

        //q now contains a valid step in the positive direction. Now we do a Line Search to find the step size that will satisfy Armijo's condition
        double epsilon = 1.0;
        double tempLogL = currentLogL;
        std::vector<double> newHyperparams;
        while (true) {
            newHyperparams = currentValue.params;
            for(int i = 0; i < numParams; i++){
                newHyperparams[i] += epsilon * q[i];
            }
            hyperparams = newHyperparams;
            updateCholesky();
            tempLogL = marginalLogLikelihood();

            double condition = 0.0;
            for(int i = 0; i < numParams; i++){
                condition += grads.back()[i] * q[i];
            }

            condition *= sufficientIncrease * epsilon;
            condition += currentLogL;

            if (tempLogL >= condition) {
                break;  // Armijo condition met
            }
            
            epsilon *= 0.5;  // Shrink step size

            if (epsilon < 1e-8) {
                std::cout << "Warning: Line search failed to find sufficient increase!" << std::endl;
                break;
            }
        }

        currentValue.score = tempLogL;
        currentValue.params = newHyperparams;
        
        history.push_back(newHyperparams);
        grads.push_back(marginalLogLikelihoodGradient());

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
        #if LOGGING==1
        if(numIter % 1 == 0){
            std::cout << "L-BFGS Optimization " << numIter << ": " << tempLogL << std::endl;
        }
        #endif
        if(numIter > 10000){
            break;
            std::cout << "Warning: L-BFGS failed to converge!" << std::endl;
        }
        numIter++;
    }
    while(converged == false);


    #if LOGGING==1
    std::cout << "Optimization has converged with parameters:" << std::endl;
    for(int i = 0; i < numParams; i++){
        std::cout << "\t" << i << ": " << std::abs(hyperparams[i]) << std::endl;
    }
    #endif

}

double BayesianOptimizer::smoothAverageAutocorrelation(const std::vector<std::vector<double>>& r){
    int L = r.size();

    assert(L >= 25);

    double sum = 0;
    for(int i = 25; i < L; i++){
        sum += autocorrelationScore(r, i-25, i);
    }

    return sum / (L - 25 + 1);
}