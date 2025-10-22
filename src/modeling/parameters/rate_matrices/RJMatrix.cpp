#include "RJMatrix.hpp"
#include "MatrixHelper.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

RJMatrix::RJMatrix(Settings& settings) : 
                                   currentQMatrix(305, 305, 0.0), oldQMatrix(305, 305, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), rLambda(settings.rLambda), omegaLambda(settings.omegaLambda), rDelta(0.5),  
                                   stationaryAlpha(30000), kDelta(0.5), omegaDelta(0.5), stationaryPriorAlpha(61, 2.0) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    currentParams[0] = Probability::Exponential::rv(&rng, kLambda);
    currentParamPriors[0] = Probability::Exponential::lnPdf(kLambda, currentParams[0]);
    currentParams[1] = Probability::Exponential::rv(&rng, rLambda);
    currentParamPriors[1] = Probability::Exponential::lnPdf(rLambda, currentParams[1]);

    for(int i = 2; i < 7; i++){
        currentParams[i] = Probability::Exponential::rv(&rng, omegaLambda);
        currentParamPriors[i] = Probability::Exponential::lnPdf(omegaLambda, currentParams[i]);
    }

    for(int i = 0; i < 7; i++){
        oldParams[i] = currentParams[i];
        oldParamPriors[i] = currentParamPriors[i];
    }

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

void RJMatrix::rebuildQMatrix(){

    double omegas[5] = {currentParams[2],0,0,0,0};
    for(int i = 1; i < 5; i++){
        omegas[i] = omegas[i-1] + currentParams[i+2];
    }

    currentQMatrix = Matrix<double>(61*currentActiveOmegas, 61*currentActiveOmegas, 0.0);
    for(auto coord : MatrixHelper::validPairs){
        for(int i = 0; i < currentActiveOmegas; i++){
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
        }
    }
    for(auto coord : MatrixHelper::nonsynonymousPairs){
        for(int i = 0; i < currentActiveOmegas; i++){   
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= omegas[i];
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= omegas[i]; 
        }
    }
    for(auto coord : MatrixHelper::transitionPairs){
        for(int i = 0; i < currentActiveOmegas; i++){
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentParams[0];
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentParams[0];
        }
    }
}

void RJMatrix::accept() {
    for(int i = 0; i < 7; i++){
        oldParams[i] = currentParams[i];
        oldParamPriors[i] = currentParamPriors[i];
    }

    oldStationary = currentStationary;
    oldStationaryPrior = currentStationaryPrior;
    oldActiveOmegas = currentActiveOmegas;

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
    else if(moveChoice == 3){
        rAcceptCount += 1;
    }

    moveChoice = -1;
}

void RJMatrix::reject() {
    for(int i = 0; i < 7; i++){
        currentParams[i] = oldParams[i];
        currentParamPriors[i] = oldParamPriors[i];
    }

    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;
    currentActiveOmegas = oldActiveOmegas;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double RJMatrix::lnPrior() {
    double prior = currentParamPriors[0] + currentStationaryPrior;

    // Include the R prior if we are swapping rates
    if(currentActiveOmegas > 1)
        prior += currentParamPriors[1];

    for(int i = 0; i < currentActiveOmegas; i++){
        prior += currentParamPriors[i + 2];
    }

    return prior;
}

double RJMatrix::updateK() {
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

double RJMatrix::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = 2;
    omegaCount += 1;

    int randomOmega = (int)(currentActiveOmegas * rng.uniformRv());

    double currentV = currentParams[2 + randomOmega];
    double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentParams[2 + randomOmega] =  newV;
    hastings = std::log(scale);

    currentParamPriors[2 + randomOmega] = Probability::Exponential::lnPdf(omegaLambda, currentParams[2 + randomOmega]);
    
    rebuildQMatrix();

    return hastings;
}

double RJMatrix::updateR() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = 3;
    rCount += 1;

    double currentV = currentParams[1];
    double scale = std::exp(rDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentParams[1] = newV;
    hastings = std::log(scale);

    currentParamPriors[1] = Probability::Exponential::lnPdf(rLambda, currentParams[1]);

    return hastings;
}

/**
 * @brief 
 */
double RJMatrix::updateActiveOmegas() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double splitAlpha = 5.0;
    double hastings = 0.0;
    this->dirty();

    double total = MatrixHelper::possibleSplit5[currentActiveOmegas-1] + MatrixHelper::possibleMerge5[currentActiveOmegas - 1];
    double splitProbs = (double)(MatrixHelper::possibleSplit5[currentActiveOmegas-1]) / total;
    double randomMove = rng.uniformRv() * total;

    if(randomMove < splitProbs){ // Perform a split
        double forwardProb = -std::log(total);
        double nextTotal = MatrixHelper::possibleSplit5[currentActiveOmegas] + MatrixHelper::possibleMerge5[currentActiveOmegas];
        double reverseProb = -std::log(nextTotal);
        hastings += reverseProb - forwardProb;

        int randomSplit = (int)(rng.uniformRv() * currentActiveOmegas);

        double u = Probability::Beta::rv(&rng, splitAlpha, splitAlpha);
        
        for(int j = currentActiveOmegas; j > randomSplit + 1; j--){
            currentParams[j + 2] = currentParams[j + 1];
        }

        double tempO = currentParams[randomSplit + 2];

        currentParams[randomSplit + 2] = tempO * u;
        currentParams[randomSplit + 3] = tempO * (1.0 - u);

        hastings += std::log(tempO);
        hastings -= Probability::Beta::lnPdf(splitAlpha, splitAlpha, u);

        currentActiveOmegas++;

        if(currentActiveOmegas == 2){         
            double newR = Probability::Exponential::rv(&rng, rLambda);
            currentParamPriors[1] = Probability::Exponential::lnPdf(rLambda, newR);
            hastings -= currentParamPriors[1];

            currentParams[1] = newR;
        }
    }
    else{ // Perform a merge
        double forwardProb = -std::log(total);
        double nextTotal = MatrixHelper::possibleSplit5[currentActiveOmegas - 2] + MatrixHelper::possibleMerge5[currentActiveOmegas - 2];
        double reverseProb = -std::log(nextTotal);
        hastings += reverseProb - forwardProb;

        int randomMerge = (int)(rng.uniformRv() * (currentActiveOmegas - 1));

        double o1 = currentParams[randomMerge + 2];
        double o2 = currentParams[randomMerge + 3];
        double total = o1 + o2;
        double u = o1 / total;

        double sum = currentParams[randomMerge + 2] + currentParams[randomMerge + 3];
        currentParams[randomMerge + 2] = sum;

        // We need to shift elements down depending on the random merge
        for(int j = randomMerge + 1; j < currentActiveOmegas-1; j++){
            currentParams[j + 2] = currentParams[j + 3];
        }

        hastings += Probability::Beta::lnPdf(splitAlpha, splitAlpha, u);
        hastings -= std::log(sum);

        currentActiveOmegas--;

        if(currentActiveOmegas == 1){
            hastings += Probability::Exponential::lnPdf(rLambda, currentParams[1]);
        }
    }

    rebuildQMatrix();

    return hastings;
}



double RJMatrix::updateStationary() {
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

Matrix<double> RJMatrix::Q() {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    int stateSpace = 61 * currentActiveOmegas;

    double scaler= 0.0;
    for(int i = 0; i < stateSpace; i++){
        double total = 0.0;
        for(int j = 0; j < stateSpace; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1.0;
        scaler += returnMatrix(i, i) * currentStationary[i % 61];
    }

    scaler = -1.0 / scaler;
    for (int i = 0; i < stateSpace; i++)
        for (int j = 0; j < stateSpace; j++)
            returnMatrix(i, j) *= scaler;
    
    if(currentActiveOmegas == 5){
        for(int i = 0; i < 61; i++){
            int offsets[5] = {0, 61, 122, 183, 244};
            for(int a = 0; a < 5; a++){
                for(int b = a + 1; b < 5; b++){
                    returnMatrix(i + offsets[a], i + offsets[b]) = currentParams[1];
                    returnMatrix(i + offsets[b], i + offsets[a]) = currentParams[1];
                }
                returnMatrix(i + offsets[a], i + offsets[a]) -= 4.0*currentParams[1];
            }
        }
    }
    else if(currentActiveOmegas == 4){
        for(int i = 0; i < 61; i++){
            int offsets[4] = {0, 61, 122, 183};
            for(int a = 0; a < 4; a++){
                for(int b = a + 1; b < 4; b++){
                    returnMatrix(i + offsets[a], i + offsets[b]) = currentParams[1];
                    returnMatrix(i + offsets[b], i + offsets[a]) = currentParams[1];
                }
                returnMatrix(i + offsets[a], i + offsets[a]) -= 3.0*currentParams[1];
            }
        }
    }
    else if(currentActiveOmegas == 3){
        for(int i = 0; i < 61; i++){
            int offsets[3] = {0, 61, 122};
            for(int a = 0; a < 3; a++){
                for(int b = a + 1; b < 3; b++){
                    returnMatrix(i + offsets[a], i + offsets[b]) = currentParams[1];
                    returnMatrix(i + offsets[b], i + offsets[a]) = currentParams[1];
                }
                returnMatrix(i + offsets[a], i + offsets[a]) -= 2.0*currentParams[1];
            }
        }
    }
    else if(currentActiveOmegas == 2){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = currentParams[1];
            returnMatrix(i + 61, i) = currentParams[1];
            returnMatrix(i, i) -= currentParams[1];
            returnMatrix(i + 61, i + 61) -= currentParams[1];
        }
    }

    return returnMatrix;
}

std::vector<double> RJMatrix::getStationary(){
    std::vector<double> returnStationary;

    double factor = 1.0/(double)(currentActiveOmegas);

    for(int i = 0; i < currentActiveOmegas; i++){
        for(double v : currentStationary){
            returnStationary.push_back(v * factor);
        }
    }
    
    return returnStationary;
}

std::array<double, 5> RJMatrix::dNdS(){
    Matrix<double> tempMatrix(305, 305, 0.0);
    std::array<double, 5> dNdS = {0,0,0,0,0}; // We start out with dN and then divide by dS
    std::array<double, 5> dS = {0,0,0,0,0};

    double omegas[5] = {currentParams[2],0,0,0,0};
    for(int i = 1; i < 5; i++){
        omegas[i] = omegas[i-1] + currentParams[i+1];
    }

    for(auto coord : MatrixHelper::validPairs){
        for(int i = 0; i < 5; i++){
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) = 1.0;
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) = 1.0;
        }
    }
    for(auto coord : MatrixHelper::nonsynonymousPairs){
        for(int i = 0; i < 5; i++){   
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) *= omegas[i];
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) *= omegas[i]; 
        }
    }
    for(auto coord : MatrixHelper::transitionPairs){
        for(int i = 0; i < 5; i++){
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentParams[0];
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentParams[0];
        }
    }

    for(auto coord : MatrixHelper::nonsynonymousPairs){
        for(int i = 0; i < 5; i++){
            dNdS[i] += tempMatrix(coord.first + (i*61), coord.second + (i*61)) * currentStationary[coord.first];
            dNdS[i] += tempMatrix(coord.second + (i*61), coord.first + (i*61)) * currentStationary[coord.second];
        }
    }

    for(auto coord : MatrixHelper::synonymousPairs){
        for(int i = 0; i < 5; i++){
            dS[i] += tempMatrix(coord.first + (i*61), coord.second + (i*61)) * currentStationary[coord.first];
            dS[i] += tempMatrix(coord.second + (i*61), coord.first + (i*61)) * currentStationary[coord.second];
        }
    }

    for(int i = 0; i < 5; i++){
        dNdS[i] /= dS[i];
    }

    return dNdS;
}

double RJMatrix::getOmegaRate(){
    return (double)omegaAcceptCount/(double)omegaCount;
}

double RJMatrix::getStationaryRate(){
    return (double)stationaryAcceptCount/(double)stationaryCount;
}

double RJMatrix::getRRate(){
    return (double)rAcceptCount/(double)rCount;
}

double RJMatrix::getKRate(){
    return (double)kAcceptCount/(double)kCount;
}

void RJMatrix::tune(){
    if(kCount > 0){
        double kRate = getKRate();

        if ( kRate > 0.33 ) {
            kDelta *= (1.0 + ((kRate-0.33)/0.67));
        }
        else {
            kDelta /= (2.0 - kRate/0.33);
        }
        kAcceptCount = 0;
        kCount = 0;
    }

    if(rCount > 0){
        double rRate = getRRate();

        if ( rRate > 0.33 ) {
            rDelta *= (1.0 + ((rRate-0.33)/0.67));
        }
        else {
            rDelta /= (2.0 - rRate/0.33);
        }
        rAcceptCount = 0;
        rCount = 0;
    }

    if(stationaryCount > 0){
        double stationaryRate = getStationaryRate();

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
        double omegaRate = getOmegaRate();

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