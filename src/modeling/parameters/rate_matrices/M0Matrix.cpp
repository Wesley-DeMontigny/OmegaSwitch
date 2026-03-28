#include "misc/Probability.hpp"
#include "misc/RandomVariable.hpp"
#include "misc/Settings.hpp"
#include "M0Matrix.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

M0Matrix::M0Matrix(Settings& settings) : 
                                   currentQMatrix(61, 61, 0.0), oldQMatrix(61, 61, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), omegaLambda(settings.omegaLambda), stationaryPriorAlpha(61, 2.0), randomStates(61, 0.0) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    currentParams[0] = Probability::Exponential::rv(&rng, kLambda);
    oldParams[0] = currentParams[0];
    currentParamPriors[0] = Probability::Exponential::lnPdf(kLambda, currentParams[0]);
    oldParamPriors[0] = currentParamPriors[0];

    currentParams[1] = Probability::Exponential::rv(&rng, omegaLambda);
    oldParams[1] = currentParams[1];
    currentParamPriors[1] = Probability::Exponential::lnPdf(omegaLambda, currentParams[1]);
    oldParamPriors[1] = currentParamPriors[1];
    
    Probability::Dirichlet::rv(&rng, stationaryPriorAlpha, currentStationary);
    oldStationary = currentStationary;
    currentStationaryPrior = Probability::Dirichlet::lnPdf(stationaryPriorAlpha, currentStationary);
    oldStationaryPrior = currentStationaryPrior;

    rebuildQMatrix();

    oldQMatrix = currentQMatrix.copy();

    std::iota(randomStates.begin(), randomStates.end(), 0);

    dirty();
}

void M0Matrix::rebuildQMatrix() {
    for (const auto& [c1, c2] : MatrixHelper::validPairs) {
        currentQMatrix(c1, c2) = currentStationary[c2];
        currentQMatrix(c2, c1) = currentStationary[c1];
    }
    for (const auto& [c1, c2] : MatrixHelper::transitionPairs) {
        currentQMatrix(c1, c2) *= currentParams[0];
        currentQMatrix(c2, c1) *= currentParams[0];
    }
    for (const auto& [c1, c2] : MatrixHelper::nonsynonymousPairs) {
        currentQMatrix(c1, c2) *= currentParams[1];
        currentQMatrix(c2, c1) *= currentParams[1];
    }
}

void M0Matrix::accept() {
    oldParams[0] = currentParams[0];
    oldParamPriors[0] = currentParamPriors[0];
    oldParams[1] = currentParams[1];
    oldParamPriors[1] = currentParamPriors[1];
    oldStationary = currentStationary;
    oldStationaryPrior = currentStationaryPrior;

    oldQMatrix = currentQMatrix.copy();

    if(moveChoice == MatrixMoves::K_MOVE){
        kAcceptCount += 1;
        if(countTuningEvents){
            tuningState->kStats.acceptCount += 1;
        }
    }
    else if(moveChoice == MatrixMoves::STATIONARY_MOVE){
        stationaryAcceptCount += 1;
        if(countTuningEvents){
            tuningState->stationaryStats.acceptCount += 1;
        }
    }
    else if(moveChoice == MatrixMoves::OMEGA_MOVE){
        omegaAcceptCount += 1;
        if(countTuningEvents){
            tuningState->omegaStats.acceptCount += 1;
        }
    }

    moveChoice = MatrixMoves::NO_MOVE;
}

void M0Matrix::reject() {
    currentParams[0] = oldParams[0];
    currentParamPriors[0] = oldParamPriors[0];
    currentParams[1] = oldParams[1];
    currentParamPriors[1] = oldParamPriors[1];
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = MatrixMoves::NO_MOVE;
}

double M0Matrix::lnPrior() {
    return currentParamPriors[0] + currentStationaryPrior + currentParamPriors[1];
}

double M0Matrix::updateK() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = MatrixMoves::K_MOVE;
    kCount += 1;
    if(countTuningEvents){
        tuningState->kStats.count += 1;
    }

    double currentV = currentParams[0];
    double scale = std::exp(tuningState->kDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentParams[0] = newV;
    hastings = std::log(scale);

    this->dirty();

    currentParamPriors[0] = Probability::Exponential::lnPdf(kLambda, currentParams[0]);

    rebuildQMatrix();

    return hastings;
}

double M0Matrix::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = MatrixMoves::OMEGA_MOVE;
    omegaCount += 1;
    if(countTuningEvents){
        tuningState->omegaStats.count += 1;
    }

    double currentV = currentParams[1];
    double scale = std::exp(tuningState->omegaDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentParams[1] = newV;
    hastings = std::log(scale);

    this->dirty();

    currentParamPriors[1] = Probability::Exponential::lnPdf(omegaLambda, currentParams[1]);

    rebuildQMatrix();

    return hastings;
}

double M0Matrix::updateStationary() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    this->dirty();
    double hastings = 0.0;

    int numElements = 1;
    moveChoice = MatrixMoves::STATIONARY_MOVE;
    stationaryCount += 1;
    if(countTuningEvents){
        tuningState->stationaryStats.count += 1;
    }

    std::vector<int> drawSet(randomStates);
    std::vector<int> randomIndices;

    for(int c = 0; c < numElements; c++){
        int i = (int)(rng.uniformRv() * drawSet.size());
        randomIndices.push_back(drawSet[i]);
        drawSet.erase(drawSet.begin() + i);
    }

    std::vector<double> x(numElements + 1, 0.0);
    std::vector<double> alphaForward(numElements + 1, 0.0);
    std::vector<double> alphaReverse(numElements + 1, 0.0);
    std::vector<double> z(numElements + 1, 0.0);


    for(int i = 0; i < 61; i++) {
        auto it = std::find(randomIndices.begin(), randomIndices.end(), i);
        if(it != randomIndices.end()) {
            x[it - randomIndices.begin()] += currentStationary[i];
        }
        else {
            x[numElements] += currentStationary[i];
        }
    }

    for(int i = 0; i < x.size(); i++) {
        alphaForward[i] = (x[i] * tuningState->stationaryAlpha) + 1.0;
    }
    
    Probability::Dirichlet::rv(&rng, alphaForward, z);

    for(int i = 0; i < z.size(); i++) {
        alphaReverse[i] = (z[i] * tuningState->stationaryAlpha) + 1.0;
    }

    double factor = z[z.size()-1] / x[x.size()-1];
    double sum = 0.0;
    for(int i = 0; i < 61; i++) {
        auto it = std::find(randomIndices.begin(), randomIndices.end(), i);
        if(it != randomIndices.end()) {
            currentStationary[i] = z[it - randomIndices.begin()];
        }
        else {
            currentStationary[i] = currentStationary[i] * factor;
        }

        sum += currentStationary[i];
    }

    // Try to rescale to avoid things shrinking to zero
    for(int i = 0; i < 61; i++) {
        currentStationary[i] = currentStationary[i]/sum;

        if(currentStationary[i] < 1E-10) {
            return -1 * INFINITY;
        }
    }

    hastings  = Probability::Dirichlet::lnPdf(alphaReverse, x) - Probability::Dirichlet::lnPdf(alphaForward, z);
    hastings += (60 - numElements) * log(factor);

    rebuildQMatrix();

    currentStationaryPrior = Probability::Dirichlet::lnPdf(stationaryPriorAlpha, currentStationary);

    return hastings;
}

Matrix<double> M0Matrix::Q() const {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    double scaler= 0.0;
    for(int i = 0; i < 61; i++){
        double total = 0.0;
        for(int j = 0; j < 61; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1;
        scaler += returnMatrix(i, i) * currentStationary[i];
    }

    scaler = -1.0 / scaler;
    for (int i = 0; i < 61; i++)
        for (int j = 0; j < 61; j++)
            returnMatrix(i, j) *= scaler;

    return returnMatrix;
}

std::vector<double> M0Matrix::getStationary(){
    return currentStationary;
}

double M0Matrix::dNdS() const {
        Matrix<double> tempMatrix(currentQMatrix.copy());

    for(const auto& [c1, c2] : MatrixHelper::validPairs){
        tempMatrix(c1, c2) /= currentStationary[c2];
        tempMatrix(c2, c1) /= currentStationary[c1];
    }

    double dN1 = 0.0;
    for(const auto& [c1, c2] : MatrixHelper::nonsynonymousPairs){
        dN1 += tempMatrix(c1, c2) * currentStationary[c1];
        dN1 += tempMatrix(c2, c1) * currentStationary[c2];
    }

    double dS1 = 0.0;
    for(const auto& [c1, c2] : MatrixHelper::synonymousPairs){
        dS1 += tempMatrix(c1, c2) * currentStationary[c1];
        dS1 += tempMatrix(c2, c1) * currentStationary[c2];
    }

    return dN1/dS1;
}

double M0Matrix::getOmegaRate() const {
    return (double)omegaAcceptCount/(double)omegaCount;
}

double M0Matrix::getStationaryRate() const {
    return (double)stationaryAcceptCount/(double)stationaryCount;
}

double M0Matrix::getKRate() const {
    return (double)kAcceptCount/(double)kCount;
}

void M0Matrix::tune(){
    if(tuningState->kStats.count > 0){
        double kRate = (double)tuningState->kStats.acceptCount/(double)tuningState->kStats.count;

        if(kRate > 0.33){
            tuningState->kDelta *= (1.0 + ((kRate-0.33)/0.67));
        }
        else{
            tuningState->kDelta /= (2.0 - kRate/0.33);
        }
        tuningState->kStats.acceptCount = 0;
        tuningState->kStats.count = 0;
    }

    if(tuningState->stationaryStats.count > 0){
        double stationaryRate = (double)tuningState->stationaryStats.acceptCount/(double)tuningState->stationaryStats.count;

        if(stationaryRate > 0.33){
            tuningState->stationaryAlpha /= (1.0 + ((stationaryRate-0.33)/0.67));
        }
        else{
            tuningState->stationaryAlpha *= (2.0 - stationaryRate/0.33);
        }

        tuningState->stationaryStats.acceptCount = 0;
        tuningState->stationaryStats.count = 0;
    }

    if(tuningState->omegaStats.count > 0){
        double omegaRate = (double)tuningState->omegaStats.acceptCount/(double)tuningState->omegaStats.count;

        if(omegaRate > 0.33){
            tuningState->omegaDelta *= (1.0 + ((omegaRate-0.33)/0.67));
        }
        else{
            tuningState->omegaDelta /= (2.0 - omegaRate/0.33);
        }
        tuningState->omegaStats.acceptCount = 0;
        tuningState->omegaStats.count = 0;
    }
}
