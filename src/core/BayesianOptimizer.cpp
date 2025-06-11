#include "BayesianOptimizer.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>
#include "Math.hpp"
#include "Matrix.hpp"
#include "RandomVariable.hpp"

BayesianOptimizer::BayesianOptimizer(int nP, int s) : numParams(nP), iterationsPerSample(s), hyperparams(nP, 1.0), defaultHyperparams(nP, 1.0) {

}

BayesianOptimizer::~BayesianOptimizer() {}

std::vector<double> BayesianOptimizer::getMaximum(){
    double max = -1.0;
    int maxIndex = -1;
    for(int i = 0; i < samples.size(); i++){
        if(objectives[i] > max){
            max = objectives[i];
            maxIndex = i;
        }
    }

    return samples[maxIndex];
}

double BayesianOptimizer::autocorrelationScore(const std::vector<double>& r, int end){
    assert(end <= r.size());

    double mean = 0.0;
    for(int i = 0; i < end; i++)
        mean += r[i];
    mean /= end;

    double variance = 0.0;
    for(int i = 0; i < end; i++)
        variance += std::pow(r[i] - mean, 2);
    variance /= end;

    int lMax = end;
    double autoCorrSum = 0.0;

    for(int i = 1; i < lMax; i++){
        double sum = 0;
        for(int t = 0; t < end - i; t++){
            sum += (r[t] - mean) * (r[t+i] - mean);
        }
        sum /= (end - i)*variance;

        autoCorrSum += std::abs(sum);
    }

    return 1 - (autoCorrSum / lMax - 1);
}

double BayesianOptimizer::kernel(std::vector<double>& a, std::vector<double>& b){
    double sum = 0.0;
    for(int i = 0; i < numParams; i++){
        double diff = a[i] - b[i];
        sum += std::pow(diff, 2) / std::pow(hyperparams[i], 2);
    }
    sum /= -2.0;

    return std::exp(sum);
}

double BayesianOptimizer::marginalLikelihood(){
    int numSamples = samples.size();
    double likelihood = -0.5;

    double sum = 0.0;
    double logDet = 0.0;
    for(int i = 0; i < numSamples; i++){
        sum += objectives[i] * alpha[i];
        logDet += std::log(choleskyFactor(i,i));
    }
    likelihood *= sum;
    likelihood -= logDet;

    likelihood -= numSamples * 0.5 * std::log(M_PI * 2);

    return likelihood;
}

std::vector<double> BayesianOptimizer::maximizeAcquisition(){
    assert(samples.size() >= 3); // We need to make sure that we have selected three coherent points

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    int numLHCSamples = 300;
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

            lhcSamples[assigningIndex][d] = increment*n + (rng.uniformRv() * increment);
        }
    }

    int maxIndex = 0;
    double maxExpectedImprovement = -1;
    for(int i = 0; i < numLHCSamples; i++){
        std::vector<double> kStar(numSamples, 0.0);
        for(int j = 0; j < numSamples; j++){
            kStar[j] = kernel(lhcSamples[i], samples[j]);
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

        // Compute expected improvement and compare
    }
}

double BayesianOptimizer::objective(std::vector<double> r) {
    int L = r.size();

    assert(L >= 25);

    double sum = 0;
    for(int i = 25; i < L; i++){
        sum += autocorrelationScore(r, i);
    }

    return sum / (L - 25 + 1);
}

void BayesianOptimizer::setBounds(std::vector<double>& sd, std::vector<double>& m){
    assert(sd.size() == m.size());
    upperBound.clear();
    lowerBound.clear();

    for(int i = 0; i < sd.size(); i++){
        upperBound.push_back(m[i] + 2.5*sd[i]);
        lowerBound.push_back(std::max(m[i] - 2.5*sd[i], 1e-3));
        defaultHyperparams[i] = (m[i] - upperBound[i]) * 0.1;
    }
}

// See Rasmussen and Williams (2008) Algorithm 2.1
void BayesianOptimizer::updateCholesky(){
    int numSamples = samples.size();

    Matrix<double> K(numSamples, numSamples, 0.0);
    for(int i = 0; i < numSamples; i++){
        for(int j = 0; j < numSamples; j++){
            K(i, j) = kernel(samples[i], samples[j]);
            if(i == j)
                K(i, j) += 1e-6; // For numerical stability
        }
    }

    choleskyFactor = Matrix<double>(numSamples, numSamples, 0.0);
    Math::choleskyDecomposition(K, choleskyFactor);
    
    alpha = objectives; // This will be overwritten
    Math::forwardSubstitutionRow(choleskyFactor, alpha);

    Matrix<double> LT;
    Math::transposeMatrix(choleskyFactor, LT);
    Math::backSubstitutionRow(LT, alpha);
}

void BayesianOptimizer::updateGaussianProcess() {
    assert(samples.size() >= 3); // We need to make sure that we have selected three coherent points

    hyperparams = defaultHyperparams; // Default hyperparams based purely on the standard deviation in acceptance rate updating
    // TODO: Optimize hyperparameters here
    updateCholesky();
}