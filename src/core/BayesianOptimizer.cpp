#include "BayesianOptimizer.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>
#include "Math.hpp"
#include "Matrix.hpp"
#include "RandomVariable.hpp"
#include "Probability.hpp"


struct SimplexVertex {
    std::vector<double> params;
    double logLikelihood;

    bool operator<(const SimplexVertex& other) const {
        return logLikelihood < other.logLikelihood;
    }
};

BayesianOptimizer::BayesianOptimizer(int nP, int s) : numParams(nP), iterationsPerSample(s), hyperparams(nP, 1.0) {

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

double BayesianOptimizer::marginalLogLikelihood(){
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

    double currentMaxObjective = 0.0;
    for(int i = 0; i < numSamples; i++){
        if(objectives[i] > currentMaxObjective)
            currentMaxObjective = objectives[i];
    }

    #if LOGGING==1
    std::cout << "Maximizing expected improvement..." << std::endl;
    #endif

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
        double predictive_std = std::sqrt(predictive_variance);

        // Compute expected improvement and compare
        double diff = predictive_mean - currentMaxObjective;
        double z = diff/predictive_std;
        double EI = std::max(0.0, diff) + 
                    predictive_std * Probability::Normal::pdf(0, 1, z) - 
                    std::abs(diff) * Probability::Normal::cdf(0, 1, z);
        
        if(maxExpectedImprovement < EI){
            maxExpectedImprovement = EI;
            maxIndex = i;
        }
    }

    #if LOGGING==1
    std::cout << "Expected improvement maximized at " << maxExpectedImprovement << std::endl;
    #endif

    return lhcSamples[maxIndex];
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

void BayesianOptimizer::setBounds(std::vector<double>& diff, std::vector<double>& mean, std::vector<bool>& sign){
    assert(diff.size() == mean.size() && sign.size() == mean.size());
    upperBound.clear();
    lowerBound.clear();

    for(int i = 0; i < diff.size(); i++){
        if(sign[i]){
            upperBound.push_back(mean[i] + diff[i]);
            lowerBound.push_back(std::max(mean[i] - 3*diff[i], 1e-3));
            hyperparams[i] = std::abs(mean[i] - upperBound[i]) * 0.1;
        }
        else {
            upperBound.push_back(mean[i] + 3*diff[i]);
            lowerBound.push_back(std::max(mean[i] - diff[i], 1e-3));
            hyperparams[i] = std::abs(mean[i] - lowerBound[i]) * 0.1;
        }
    }

    #if LOGGING==1
    std::cout << "Setting bounds on the hyperparameter optimization:" << std::endl;
    for(int i = 0; i < upperBound.size(); i++){
        std::cout << "\t" << i << ": " << lowerBound[i] << "-" << upperBound[i] << std::endl;
    }
    std::cout << "Setting intial hyperparameters:" << std::endl;
    for(int i = 0; i < hyperparams.size(); i++){
        std::cout << "\t" << i << ": " << hyperparams[i] << std::endl;
    }
    #endif
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

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // Initialize simplex
    updateCholesky();
    std::vector<SimplexVertex> hyperparamPoints;
    SimplexVertex initPoint;
    initPoint.params = hyperparams;
    initPoint.logLikelihood = marginalLogLikelihood();
    hyperparamPoints.push_back(initPoint);

    for(int i = 0; i < numParams; i++){
        SimplexVertex newPoint;
        newPoint.params = hyperparamPoints[0].params;
        for(int j = 0; j < numParams; j++){
            double modifiedDim = newPoint.params[j] * std::exp(0.25*(rng.uniformRv() - 0.5));
            modifiedDim = std::clamp(modifiedDim, lowerBound[j], upperBound[j]);
            newPoint.params[j] = modifiedDim;
        }
        hyperparams = newPoint.params;
        updateCholesky();
        newPoint.logLikelihood = marginalLogLikelihood();
        hyperparamPoints.push_back(newPoint);
    }

    std::sort(hyperparamPoints.begin(), hyperparamPoints.end(), [](const SimplexVertex& a, const SimplexVertex& b) {
        return a.logLikelihood > b.logLikelihood;
    });

    #if LOGGING==1
    std::cout << "Initialized hyperparameter simplex..." << std::endl;
    #endif


    double reflectParam = 1.0;
    double expandParam = 2.0;
    double contractParam = 0.5;
    double shrinkParam = 0.5;

    int numSimplexPoints = hyperparamPoints.size();
    int iterations = 0;
    int maxIterations = 1000;

    #if LOGGING==1
    std::cout << "Starting Nelder-Mead optimization..." << std::endl;
    #endif

    bool converged = false;
    do{
        if(iterations >= maxIterations){
            std::cout << "Warning: Nelder-Mead Optimization Failed to Converge in " << maxIterations << " Iterations..." << std::endl;
            break;
        }
        iterations++;
        // Compute the centroid
        std::vector<double> centroid(numParams, 0.0);
        for(int i = 1; i < numSimplexPoints; i++){
            for(int j = 0; j < numParams; j++){
                centroid[j] += hyperparamPoints[i].params[j];
            }
        }
        for(int i = 0; i < numParams; i++){
            centroid[i] /= numSimplexPoints - 1;
        }

        // Do we reflect the point?
        std::vector<double> reflected = centroid;
        for(int i = 0; i < numParams; i++){
            reflected[i] += reflectParam * (centroid[i] - hyperparamPoints[numSimplexPoints-1].params[i]);
        }
        hyperparams = reflected;
        updateCholesky();
        double reflectedMarginalLL = marginalLogLikelihood();
        if(reflectedMarginalLL > hyperparamPoints[numSimplexPoints-1].logLikelihood){
            // Can we get away with an expansion?
            if(reflectedMarginalLL > hyperparamPoints[0].logLikelihood){
                std::vector<double> expanded = centroid;
                for(int i = 0; i < numParams; i++){
                    expanded[i] += expandParam * (reflected[i] - centroid[i]);
                }
                hyperparams = expanded;
                updateCholesky();
                double expandedMarginalLL = marginalLogLikelihood();
                if(expandedMarginalLL > reflectedMarginalLL){
                    // Accept expansion
                    hyperparamPoints[numSimplexPoints-1].logLikelihood = expandedMarginalLL;
                    hyperparamPoints[numSimplexPoints-1].params = expanded;
                }
                else {
                    // Accept reflection
                    hyperparamPoints[numSimplexPoints-1].logLikelihood = reflectedMarginalLL;
                    hyperparamPoints[numSimplexPoints-1].params = reflected;
                }
            }
            else{
                // Accept reflection
                hyperparamPoints[numSimplexPoints-1].logLikelihood = reflectedMarginalLL;
                hyperparamPoints[numSimplexPoints-1].params = reflected;
            }
        }
        else {
            // Should we contract the point?
            std::vector<double> contracted  = centroid;
            for(int i = 0; i < numParams; i++){
                contracted[i] += contractParam * (reflected[i] - centroid[i]);
            }
            hyperparams = contracted;
            updateCholesky();
            double contractedMarginalLL = marginalLogLikelihood();
            if(contractedMarginalLL > reflectedMarginalLL){
                // Accept contraction
                hyperparamPoints[numSimplexPoints-1].logLikelihood = contractedMarginalLL;
                hyperparamPoints[numSimplexPoints-1].params = contracted;
            }
            else{ // Shrink the whole simplex towards the max
                for(int i = 1; i < numSimplexPoints; i++){
                    for(int j = 0; j < numParams; j++){
                        hyperparamPoints[i].params[j] = hyperparamPoints[0].params[j] - shrinkParam * (hyperparamPoints[i].params[j] - hyperparamPoints[0].params[j]);
                    }
                    hyperparams = hyperparamPoints[i].params;
                    updateCholesky();
                    hyperparamPoints[i].logLikelihood = marginalLogLikelihood();
                }
            }
        }

        std::sort(hyperparamPoints.begin(), hyperparamPoints.end(), [](const SimplexVertex& a, const SimplexVertex& b) {
            return a.logLikelihood > b.logLikelihood;
        });

        double max_diff = 0.0;
        for(int i = 1; i < hyperparamPoints.size(); i++){
            for(int j = 0; j < numParams; j++){
                double diff = abs(hyperparamPoints[i].params[j] - hyperparamPoints[0].params[j]);
                if(diff > max_diff){
                    max_diff = diff;
                }
            }
        }

        if(max_diff < 1e-2)
            converged = true;
    }
    while(converged == false);


    hyperparams = hyperparamPoints[0].params;
    updateCholesky();
    
    #if LOGGING==1
    std::cout << "Optimization has convered with parameters:" << std::endl;
    for(int i = 0; i < numParams; i++){
        std::cout << "\t" << i << ": " << hyperparams[i] << std::endl;
    }
    #endif
}