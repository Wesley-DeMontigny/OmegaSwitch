#include "M0Matrix.hpp"
#include "MatrixHelper.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

M0Matrix::M0Matrix(Settings& settings) : 
                                   currentQMatrix(61, 61, 0.0), oldQMatrix(61, 61, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), omegaLambda(settings.omegaLambda), stationaryAlpha(30000), kDelta(0.5),
                                   omegaDelta(0.5), stationaryPriorAlpha(61, 2.0) {

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

    for(int i = 0; i < 61; i++)
        randomStates.push_back(i);

    dirty();
}

void M0Matrix::rebuildQMatrix() {
    for (auto coord : MatrixHelper::validPairs) {
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.second];
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.first];
    }
    for (auto coord : MatrixHelper::transitionPairs) {
        currentQMatrix(coord.first, coord.second) *= currentParams[0];
        currentQMatrix(coord.second, coord.first) *= currentParams[0];
    }
    for (auto coord : MatrixHelper::nonsynonymousPairs) {
        currentQMatrix(coord.first, coord.second) *= currentParams[1];
        currentQMatrix(coord.second, coord.first) *= currentParams[1];
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

    if(moveChoice == 0){
        kAcceptCount += 1;
    }
    else if(moveChoice == 1){
        stationaryAcceptCount += 1;
    }
    else if(moveChoice == 2){
        omegaAcceptCount += 1;
    }

    moveChoice = -1;
}

void M0Matrix::reject() {
    currentParams[0] = oldParams[0];
    currentParamPriors[0] = oldParamPriors[0];
    currentParams[1] = oldParams[1];
    currentParamPriors[1] = oldParamPriors[1];
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double M0Matrix::lnPrior() {
    return currentParamPriors[0] + currentStationaryPrior + currentParamPriors[1];
}

double M0Matrix::updateK() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 0;
    kCount += 1;

    double currentV = currentParams[0];
    double scale = std::exp(kDelta * (rng.uniformRv() - 0.5));
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

    moveChoice = 2;
    omegaCount += 1;

    double currentV = currentParams[1];
    double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
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

    int numElements = 30;
    moveChoice = 1;
    stationaryCount += 1;

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
        alphaForward[i] = (x[i] * stationaryAlpha) + 1.0;
    }
    
    Probability::Dirichlet::rv(&rng, alphaForward, z);

    for(int i = 0; i < z.size(); i++) {
        alphaReverse[i] = (z[i] * stationaryAlpha) + 1.0;
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

Matrix<double> M0Matrix::Q() {
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

double M0Matrix::dNdS(){
        Matrix<double> tempMatrix(currentQMatrix.copy());

    for(auto coord : MatrixHelper::validPairs){
        tempMatrix(coord.first, coord.second) /= currentStationary[coord.second];
        tempMatrix(coord.second, coord.first) /= currentStationary[coord.first];
    }

    double dN1 = 0.0;
    for(auto coord : MatrixHelper::nonsynonymousPairs){
        dN1 += tempMatrix(coord.first, coord.second) * currentStationary[coord.first];
        dN1 += tempMatrix(coord.second, coord.first) * currentStationary[coord.second];
    }

    double dS1 = 0.0;
    for(auto coord : MatrixHelper::synonymousPairs){
        dS1 += tempMatrix(coord.first, coord.second) * currentStationary[coord.first];
        dS1 += tempMatrix(coord.second, coord.first) * currentStationary[coord.second];
    }

    return dN1/dS1;
}

double M0Matrix::getOmegaRate(){
    return (double)omegaAcceptCount/(double)omegaCount;
}

double M0Matrix::getStationaryRate(){
    return (double)stationaryAcceptCount/(double)stationaryCount;
}

double M0Matrix::getKRate(){
    return (double)kAcceptCount/(double)kCount;
}

void M0Matrix::tune(){
    if(kCount > 0){
        double kRate = (double)kAcceptCount/(double)kCount;

        if ( kRate > 0.33 ) {
            kDelta *= (1.0 + ((kRate-0.33)/0.67));
        }
        else {
            kDelta /= (2.0 - kRate/0.33);
        }
        kAcceptCount = 0;
        kCount = 0;
    }

    if(stationaryCount > 0){
        double stationaryRate = (double)stationaryAcceptCount/(double)stationaryCount;

        if ( stationaryRate > 0.33 ) {
            stationaryAlpha /= (1.0 + ((stationaryRate-0.33)/0.67));
        }
        else {
            stationaryAlpha *= (2.0 - stationaryRate/0.33);
        }

        stationaryAcceptCount = 0;
        stationaryCount = 0;
    }

    if(omegaCount > 0){
        double omegaRate = (double)omegaAcceptCount/(double)omegaCount;

        if ( omegaRate > 0.33 ) {
            omegaDelta *= (1.0 + ((omegaRate-0.33)/0.67));
        }
        else {
            omegaDelta /= (2.0 - omegaRate/0.33);
        }
        omegaAcceptCount = 0;
        omegaCount = 0;
    }
}