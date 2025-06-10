#include "BayesianOptimizer.hpp"
#include <cmath>
#include <cassert>
#include "Math.hpp"
#include "Matrix.hpp"
#include "EigenSystem.hpp"

BayesianOptimizer::BayesianOptimizer(int nP, int s) : numParams(nP), iterationsPerSample(s), hyperparams(nP, 1.0) {
}

BayesianOptimizer::~BayesianOptimizer() {}

double BayesianOptimizer::objective(std::vector<double> r) {
    int L = r.size();

    assert(L >= 25);

    double sum = 0;
    for(int i = 25; i < L; i++){
        sum += autocorrelationScore(r, i);
    }

    return sum / (L - 25 + 1);
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

double BayesianOptimizer::kernel(std::vector<double> a, std::vector<double> b){
    double sum = 0.0;
    for(int i = 0; i < numParams; i++){
        double diff = a[i] - b[i];
        sum += std::pow(diff, 2) / std::pow(hyperparams[i], 2);
    }
    sum /= -2.0;

    return std::exp(sum);
}