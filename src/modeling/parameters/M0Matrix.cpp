#include "M0Matrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

M0Matrix::M0Matrix(Settings settings) : 
                                   currentQMatrix(61, 61, 0.0), oldQMatrix(61, 61, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), omegaLambda(settings.omegaLambda), stationaryAlpha(75000), kDelta(0.25),
                                   omegaDelta(0.25), stationaryPriorAlpha(61, 2.0) {
    std::vector<int> aaMap = {8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 19, 19, 15, 15, 15, 15, 1, 18, 1, 9, 4, 9, 4};  
    std::vector<const char*> codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};


    // Because of the complicated nature of this matrix, we need to classify each of the positions in the matrix;
    for(int i = 0; i < 61; i++){
        for(int j = i + 1; j < 61; j++){
            int mismatch = 0;
            bool isTransition = false;
            for(int k = 0; k < 3; k++){
                if(codons[i][k] != codons[j][k]){
                    mismatch++;
                    if(mismatch > 1){
                        break;
                    }
                    if((codons[i][k] == 'A' && codons[j][k] == 'G') || (codons[i][k] == 'G' && codons[j][k] == 'A') || 
                       (codons[i][k] == 'T' && codons[j][k] == 'C') || (codons[i][k] == 'C' && codons[j][k] == 'T'))
                        isTransition = true;
                }
            }
            if(mismatch == 1){
                auto pair = std::make_pair(i, j);
                valid.insert(pair);
                if(aaMap[i] != aaMap[j])
                    nonsynonymous.insert(pair);
                else
                    synonymous.insert(pair);
                if(isTransition)
                    transition.insert(pair);
            }
        }
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    currentK = Probability::Exponential::rv(&rng, kLambda);
    oldK = currentK;
    currentKPrior = Probability::Exponential::lnPdf(kLambda, currentK);
    oldKPrior = currentKPrior;

    currentOmega = Probability::Exponential::rv(&rng, omegaLambda);
    oldOmega = currentOmega;
    currentOmegaPrior = Probability::Exponential::lnPdf(omegaLambda, currentOmega);
    oldOmegaPrior = currentOmegaPrior;
    
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
    for (auto coord : valid) {
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.second];
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.first];
    }
    for (auto coord : transition) {
        currentQMatrix(coord.first, coord.second) *= currentK;
        currentQMatrix(coord.second, coord.first) *= currentK;
    }
    for (auto coord : nonsynonymous) {
        currentQMatrix(coord.first, coord.second) *= currentOmega;
        currentQMatrix(coord.second, coord.first) *= currentOmega;
    }
}

void M0Matrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;
    oldOmega = currentOmega;
    oldOmegaPrior = currentOmegaPrior;
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
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentOmega = oldOmega;
    currentOmegaPrior = oldOmegaPrior;
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double M0Matrix::lnPrior() {
    return currentKPrior + currentStationaryPrior + currentOmegaPrior;
}

double M0Matrix::updateK() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 0;
    kCount += 1;

    double currentV = currentK;
    double scale = std::exp(kDelta * (rng.uniformRv() - 0.5));
    double newV = currentK * scale;

    currentK = newV;
    hastings = std::log(scale);

    this->dirty();

    currentKPrior = Probability::Exponential::lnPdf(kLambda, currentK);

    rebuildQMatrix();

    return hastings;
}

double M0Matrix::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 2;
    omegaCount += 1;

    double currentV = currentOmega;
    double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
    double newV = currentOmega * scale;

    currentOmega = newV;
    hastings = std::log(scale);

    this->dirty();

    currentOmegaPrior = Probability::Exponential::lnPdf(omegaLambda, currentOmega);

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

void M0Matrix::tune(){
    double kRate = (double)kAcceptCount/(double)kCount;

    if ( kRate > 0.33 ) {
        kDelta *= (1.0 + ((kRate-0.33)/0.67));
    }
    else {
        kDelta /= (2.0 - kRate/0.33);
    }

    kAcceptCount = 0;
    kCount = 0;

    double omegaRate = (double)omegaAcceptCount/(double)omegaCount;

    if ( omegaRate > 0.33 ) {
        omegaDelta *= (1.0 + ((omegaRate-0.33)/0.67));
    }
    else {
        omegaDelta /= (2.0 - omegaRate/0.33);
    }

    omegaAcceptCount = 0;
    omegaCount = 0;

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