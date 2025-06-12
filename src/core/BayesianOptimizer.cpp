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

double BayesianOptimizer::UBC(std::vector<double>& sample, int numSamples){
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

    std::cout << "Using latin hypercube sampling to sample candidate parameters..." << std::endl;

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

    std::cout << "Maximizing acquisition function..." << std::endl;

    std::vector<ParamScorePair> parameterScorePairs;
    for(int i = 0; i < numLHCSamples; i++){
        
        //double score = expectedImprovement(lhcSamples[i], currentMaxObjective, numSamples);
        double score = UBC(lhcSamples[i], numSamples);

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
                K(i, j) += 0.1; // For numerical stability
        }
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

void BayesianOptimizer::updateGaussianProcess() {
    assert(samples.size() >= 3); // We need to make sure that we have selected three coherent points

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // Initialize simplex using latin hypercube sampling
    std::cout << "Using latin hypercube sampling to get initial Nelder-Mead simplex..." << std::endl;
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

    std::sort(hyperparamPoints.begin(), hyperparamPoints.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
        return a.score > b.score;
    });

    hyperparamPoints.erase(hyperparamPoints.begin() + numParams + 1, hyperparamPoints.end());

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

    std::cout << "Starting Nelder-Mead optimization of hyperparameters..." << std::endl;

    bool converged = false;
    do{
        if(iterations >= maxIterations){
            std::cout << "Warning: Nelder-Mead Optimization Failed to Converge in " << maxIterations << " Iterations!" << std::endl;
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
        hyperparams = pow10Vector(reflected);
        updateCholesky();
        double reflectedMarginalLL = marginalLogLikelihood();
        if(reflectedMarginalLL > hyperparamPoints[numSimplexPoints-1].score){
            // Can we get away with an expansion?
            if(reflectedMarginalLL > hyperparamPoints[0].score){
                std::vector<double> expanded = centroid;
                for(int i = 0; i < numParams; i++){
                    expanded[i] += expandParam * (reflected[i] - centroid[i]);
                }
                hyperparams = pow10Vector(expanded);
                updateCholesky();
                double expandedMarginalLL = marginalLogLikelihood();
                if(expandedMarginalLL > reflectedMarginalLL){
                    // Accept expansion
                    hyperparamPoints[numSimplexPoints-1].score = expandedMarginalLL;
                    hyperparamPoints[numSimplexPoints-1].params = expanded;
                }
                else {
                    // Accept reflection
                    hyperparamPoints[numSimplexPoints-1].score = reflectedMarginalLL;
                    hyperparamPoints[numSimplexPoints-1].params = reflected;
                }
            }
            else{
                // Accept reflection
                hyperparamPoints[numSimplexPoints-1].score = reflectedMarginalLL;
                hyperparamPoints[numSimplexPoints-1].params = reflected;
            }
        }
        else {
            // Should we contract the point?
            std::vector<double> contracted  = centroid;
            for(int i = 0; i < numParams; i++){
                contracted[i] += contractParam * (reflected[i] - centroid[i]);
            }
            hyperparams = pow10Vector(contracted);
            updateCholesky();
            double contractedMarginalLL = marginalLogLikelihood();
            if(contractedMarginalLL > reflectedMarginalLL){
                // Accept contraction
                hyperparamPoints[numSimplexPoints-1].score = contractedMarginalLL;
                hyperparamPoints[numSimplexPoints-1].params = contracted;
            }
            else{ // Shrink the whole simplex towards the max
                for(int i = 1; i < numSimplexPoints; i++){
                    for(int j = 0; j < numParams; j++){
                        hyperparamPoints[i].params[j] = hyperparamPoints[0].params[j] - shrinkParam * (hyperparamPoints[i].params[j] - hyperparamPoints[0].params[j]);
                    }
                    hyperparams = pow10Vector(hyperparamPoints[i].params);
                    updateCholesky();
                    hyperparamPoints[i].score = marginalLogLikelihood();
                }
            }
        }

        std::sort(hyperparamPoints.begin(), hyperparamPoints.end(), [](const ParamScorePair& a, const ParamScorePair& b) {
            return a.score > b.score;
        });

        double max_diff = 0.0;
        for(int i = 1; i < hyperparamPoints.size(); i++){
            for(int j = 0; j < numParams; j++){
                double diff = std::abs(hyperparamPoints[i].params[j] - hyperparamPoints[0].params[j]);
                if(diff > max_diff){
                    max_diff = diff;
                }
            }
        }

        if(iterations % 25 == 0){
            std::cout << "Nelder-Mead Optimization " << iterations << ": " << hyperparamPoints[0].score << std::endl;
        }

        if(max_diff < 1e-2)
            converged = true;
    }
    while(converged == false);


    hyperparams = pow10Vector(hyperparamPoints[0].params);
    updateCholesky();
    
    //#if LOGGING==1
    std::cout << "Optimization has converged with parameters:" << std::endl;
    for(int i = 0; i < numParams; i++){
        std::cout << "\t" << i << ": " << hyperparams[i] << std::endl;
    }
    //#endif
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