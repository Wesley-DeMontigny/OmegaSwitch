#include "M3S2Matrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

M3S2Matrix::M3S2Matrix(Settings settings) : 
                                   currentQMatrix(183, 183, 0.0), oldQMatrix(183, 183, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), gammaLambda(settings.gammaLambda), rLambda(settings.rLambda),
                                   omegaLambda(settings.omegaLambda), r1Delta(0.5), r2Delta(0.5), gammaDelta(0.5),
                                   stationaryAlpha(75000), kDelta(0.5), omega1Delta(0.5), omega2Delta(0.5), omega3Delta(0.5), 
                                   stationaryPriorAlpha(61, 2.0) {
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

    currentOmega1 = Probability::Exponential::rv(&rng, omegaLambda);
    oldOmega1 = currentOmega1;
    currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
    oldOmega1Prior = currentOmega1Prior;

    currentOmega2 = Probability::Exponential::rv(&rng, omegaLambda);
    oldOmega2 = currentOmega2;
    currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);
    oldOmega2Prior = currentOmega2Prior;

    currentOmega3 = Probability::Exponential::rv(&rng, omegaLambda);
    oldOmega3 = currentOmega3;
    currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);
    oldOmega3Prior = currentOmega3Prior;

    currentGamma = Probability::Exponential::rv(&rng, gammaLambda);
    oldGamma = currentGamma;
    currentGammaPrior = Probability::Exponential::lnPdf(gammaLambda, currentGamma);
    oldGammaPrior = currentGammaPrior;

    currentR1 = Probability::Exponential::rv(&rng, rLambda);
    oldR1 = currentR1;
    currentR1Prior = Probability::Exponential::lnPdf(rLambda, currentR1);
    oldR1Prior = currentR1Prior;

    currentR2 = Probability::Exponential::rv(&rng, rLambda);
    oldR2 = currentR2;
    currentR2Prior = Probability::Exponential::lnPdf(rLambda, currentR2);
    oldR2Prior = currentR2Prior;
    
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

void M3S2Matrix::rebuildQMatrix(){
    for(auto coord : valid){
        for(int i = 0; i < 3; i++){
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
        }
    }
    for(auto coord : transition){
        for(int i = 0; i < 3; i++){
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
        }
    }
    for(auto coord : nonsynonymous){
        currentQMatrix(coord.first, coord.second) *= currentOmega1;
        currentQMatrix(coord.second, coord.first) *= currentOmega1;

        currentQMatrix(coord.first + 61, coord.second + 61) *= currentOmega1 + currentOmega2;
        currentQMatrix(coord.second + 61, coord.first + 61) *= currentOmega1 + currentOmega2;

        currentQMatrix(coord.first + 122, coord.second + 122) *= currentOmega1 + currentOmega2 + currentOmega3;
        currentQMatrix(coord.second + 122, coord.first + 122) *= currentOmega1 + currentOmega2 + currentOmega3;
    }
}

void M3S2Matrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;
    oldOmega1 = currentOmega1;
    oldOmega1Prior = currentOmega1Prior;
    oldOmega2 = currentOmega2;
    oldOmega2Prior = currentOmega2Prior;
    oldOmega3 = currentOmega3;
    oldOmega3Prior = currentOmega3Prior;
    oldGamma = currentGamma;
    oldGammaPrior = currentGammaPrior;
    oldR1 = currentR1;
    oldR1Prior = currentR1Prior;
    oldR2 = currentR2;
    oldR2Prior = currentR2Prior;
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
        omega1AcceptCount += 1;
    }
    else if(moveChoice == 3){
        omega2AcceptCount += 1;
    }
    else if(moveChoice == 4){
        omega3AcceptCount += 1;
    }
    else if(moveChoice == 5){
        gammaAcceptCount= 1;
    }
    else if(moveChoice == 6){
        r1AcceptCount += 1;
    }
    else if(moveChoice == 7){
        r2AcceptCount += 1;
    }

    moveChoice = -1;
}

void M3S2Matrix::reject() {
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentOmega1 = oldOmega1;
    currentOmega1Prior = oldOmega1Prior;
    currentOmega2 = oldOmega2;
    currentOmega2Prior = oldOmega2Prior;
    currentOmega3 = oldOmega3;
    currentOmega3Prior = oldOmega3Prior;
    currentGamma = oldGamma;
    currentGammaPrior = oldGammaPrior;
    currentR1 = oldR1;
    currentR1Prior = oldR1Prior;
    currentR2 = oldR2;
    currentR2Prior = oldR2Prior;
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double M3S2Matrix::lnPrior() {
    return currentKPrior + currentStationaryPrior + currentOmega1Prior + currentOmega2Prior + currentOmega3Prior + currentGammaPrior + currentR1Prior + currentR2Prior;
}

double M3S2Matrix::updateK() {
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

double M3S2Matrix::updateOmega1() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 2;
    omega1Count += 1;

    double currentV = currentOmega1;
    double scale = std::exp(omega1Delta * (rng.uniformRv() - 0.5));
    double newV = currentOmega1 * scale;

    currentOmega1 = newV;
    hastings = std::log(scale);

    this->dirty();

    currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);

    rebuildQMatrix();

    return hastings;
}

double M3S2Matrix::updateOmega2() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 3;
    omega2Count += 1;

    double currentV = currentOmega2;
    double scale = std::exp(omega2Delta * (rng.uniformRv() - 0.5));
    double newV = currentOmega2 * scale;

    currentOmega2 = newV;
    hastings = std::log(scale);

    this->dirty();

    currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);

    rebuildQMatrix();

    return hastings;
}

double M3S2Matrix::updateOmega3() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 4;
    omega3Count += 1;

    double currentV = currentOmega3;
    double scale = std::exp(omega3Delta * (rng.uniformRv() - 0.5));
    double newV = currentOmega3 * scale;

    currentOmega3 = newV;
    hastings = std::log(scale);

    this->dirty();

    currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);

    rebuildQMatrix();

    return hastings;
}

double M3S2Matrix::updateGamma() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 5;
    gammaCount += 1;

    double currentV = currentGamma;
    double scale = std::exp(gammaDelta * (rng.uniformRv() - 0.5));
    double newV = currentGamma * scale;

    currentGamma = newV;
    hastings = std::log(scale);

    this->dirty();

    currentGammaPrior = Probability::Exponential::lnPdf(gammaLambda, currentGamma);

    return hastings;
}

double M3S2Matrix::updateR1() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 6;
    r1Count += 1;

    double currentV = currentR1;
    double scale = std::exp(r1Delta * (rng.uniformRv() - 0.5));
    double newV = currentR1 * scale;

    currentR1 = newV;
    hastings = std::log(scale);

    this->dirty();

    currentR1Prior = Probability::Exponential::lnPdf(rLambda, currentR1);

    return hastings;
}

double M3S2Matrix::updateR2() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 7;
    r2Count += 1;

    double currentV = currentR2;
    double scale = std::exp(r2Delta * (rng.uniformRv() - 0.5));
    double newV = currentR2 * scale;

    currentR2 = newV;
    hastings = std::log(scale);

    this->dirty();

    currentR2Prior = Probability::Exponential::lnPdf(rLambda, currentR2);

    return hastings;
}

double M3S2Matrix::updateStationary() {
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

Matrix<double> M3S2Matrix::Q() {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    double scaler= 0.0;
    for(int i = 0; i < 183; i++){
        double total = 0.0;
        for(int j = 0; j < 183; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1;
        scaler += returnMatrix(i, i) * currentStationary[i % 61];
    }

    scaler = -1.0 / scaler;
    for (int i = 0; i < 183; i++)
        for (int j = 0; j < 183; j++)
            returnMatrix(i, j) *= scaler;
    

    for(int i = 0; i < 61; i++){
        returnMatrix(i, i + 61) = currentGamma;
        returnMatrix(i + 61, i) = currentGamma;
        returnMatrix(i, i + 122) = currentGamma * currentR1;
        returnMatrix(i + 122, i) = currentGamma * currentR1;
        returnMatrix(i + 61, i + 122) = currentGamma * currentR2;
        returnMatrix(i + 122, i + 61) = currentGamma * currentR2;

        returnMatrix(i, i) -= currentGamma * (1 + currentR1);
        returnMatrix(i + 61, i + 61) -= currentGamma * (1 + currentR2);
        returnMatrix(i + 122, i + 122) -= currentGamma * (currentR1 + currentR2);
    }

    return returnMatrix;
}

std::vector<double> M3S2Matrix::getStationary(){
    std::vector<double> returnStationary;

    for(int i = 0; i < 3; i++){
        for(double v : currentStationary){
            returnStationary.push_back(v/3);
        }
    }
    
    return returnStationary;
}

std::tuple<double, double, double> M3S2Matrix::dNdS(){
        Matrix<double> tempMatrix(currentQMatrix.copy());

    for(auto coord : valid){
        tempMatrix(coord.first, coord.second) /= currentStationary[coord.second];
        tempMatrix(coord.second, coord.first) /= currentStationary[coord.first];

        tempMatrix(coord.first + 61, coord.second +  61) /= currentStationary[coord.second];
        tempMatrix(coord.second + 61, coord.first + 61) /= currentStationary[coord.first];

        tempMatrix(coord.first + 122, coord.second +  122) /= currentStationary[coord.second];
        tempMatrix(coord.second + 122, coord.first + 122) /= currentStationary[coord.first];
    }

    double dN1 = 0.0;
    double dN2 = 0.0;
    double dN3 = 0.0;
    for(auto coord : nonsynonymous){
        dN1 += tempMatrix(coord.first, coord.second) * currentStationary[coord.first];
        dN1 += tempMatrix(coord.second, coord.first) * currentStationary[coord.second];

        dN2 += tempMatrix(coord.first + 61, coord.second + 61) * currentStationary[coord.first];
        dN2 += tempMatrix(coord.second + 61, coord.first + 61) * currentStationary[coord.second];

        dN3 += tempMatrix(coord.first + 122, coord.second + 122) * currentStationary[coord.first];
        dN3 += tempMatrix(coord.second + 122, coord.first + 122) * currentStationary[coord.second];
    }

    double dS1 = 0.0;
    double dS2 = 0.0;
    double dS3 = 0.0;
    for(auto coord : synonymous){
        dS1 += tempMatrix(coord.first, coord.second) * currentStationary[coord.first];
        dS1 += tempMatrix(coord.second, coord.first) * currentStationary[coord.second];

        dS2 += tempMatrix(coord.first + 61, coord.second + 61) * currentStationary[coord.first];
        dS2 += tempMatrix(coord.second + 61, coord.first + 61) * currentStationary[coord.second];

        dS3 += tempMatrix(coord.first + 61, coord.second + 61) * currentStationary[coord.first];
        dS3 += tempMatrix(coord.second + 122, coord.first + 122) * currentStationary[coord.second];
    }

    return std::make_tuple(dN1/dS1, dN2/dS2, dN3/dS3);
}

void M3S2Matrix::tune(){
    double kRate = (double)kAcceptCount/(double)kCount;

    if ( kRate > 0.33 ) {
        kDelta *= (1.0 + ((kRate-0.33)/0.67));
    }
    else {
        kDelta /= (2.0 - kRate/0.33);
    }
    kAcceptCount = 0;
    kCount = 0;

    double gammaRate = (double)gammaAcceptCount/(double)gammaCount;

    if ( gammaRate > 0.33 ) {
        gammaDelta *= (1.0 + ((gammaRate-0.33)/0.67));
    }
    else {
        gammaDelta /= (2.0 - gammaRate/0.33);
    }
    gammaAcceptCount = 0;
    gammaCount = 0;

    double r1Rate = (double)r1AcceptCount/(double)r1Count;

    if ( r1Rate > 0.33 ) {
        r1Delta *= (1.0 + ((r1Rate-0.33)/0.67));
    }
    else {
        r1Delta /= (2.0 - r1Rate/0.33);
    }
    r1AcceptCount = 0;
    r1Count = 0;

    double r2Rate = (double)r2AcceptCount/(double)r2Count;

    if ( r2Rate > 0.33 ) {
        r2Delta *= (1.0 + ((r2Rate-0.33)/0.67));
    }
    else {
        r2Delta /= (2.0 - r2Rate/0.33);
    }
    r2AcceptCount = 0;
    r2Count = 0;

    double omega1Rate = (double)omega1AcceptCount/(double)omega1Count;

    if ( omega1Rate > 0.33 ) {
        omega1Delta *= (1.0 + ((omega1Rate-0.33)/0.67));
    }
    else {
        omega1Delta /= (2.0 - omega1Rate/0.33);
    }
    omega1AcceptCount = 0;
    omega1Count = 0;

    double omega2Rate = (double)omega1AcceptCount/(double)omega1Count;

    if ( omega2Rate > 0.33 ) {
        omega2Delta *= (1.0 + ((omega2Rate-0.33)/0.67));
    }
    else {
        omega2Delta /= (2.0 - omega2Rate/0.33);
    }
    omega2AcceptCount = 0;
    omega2Count = 0;

    double omega3Rate = (double)omega1AcceptCount/(double)omega1Count;

    if ( omega3Rate > 0.33 ) {
        omega3Delta *= (1.0 + ((omega3Rate-0.33)/0.67));
    }
    else {
        omega3Delta /= (2.0 - omega3Rate/0.33);
    }
    omega3AcceptCount = 0;
    omega3Count = 0;

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