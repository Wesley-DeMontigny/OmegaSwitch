#include "RJDPPMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

RJDPPMatrix::RJDPPMatrix(Settings settings) : 
                                   currentQMatrix(183, 183, 0.0), oldQMatrix(183, 183, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), rLambda(settings.rLambda), rDelta(0.5),  stationaryAlpha(30000), kDelta(0.5),
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

    currentR = Probability::Exponential::rv(&rng, rLambda);
    oldR = currentR;
    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);
    oldRPrior = currentRPrior;

    Probability::Dirichlet::rv(&rng, stationaryPriorAlpha, currentStationary);
    oldStationary = currentStationary;
    currentStationaryPrior = Probability::Dirichlet::lnPdf(stationaryPriorAlpha, currentStationary);
    oldStationaryPrior = currentStationaryPrior;

    refreshQBackground(3);
    
    oldQMatrix = currentQMatrix.copy();

    for(int i = 0; i < 61; i++)
        randomStates.push_back(i);

    dirty();
}

void RJDPPMatrix::refreshQBackground(int numClasses){
    if(numClasses == 3){
        currentQMatrix = Matrix<double>(183, 183, 0.0);
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
    }
    else if(numClasses == 2){
        currentQMatrix = Matrix<double>(122, 122, 0.0);
        for(auto coord : valid){
            for(int i = 0; i < 2; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
            }
        }
        for(auto coord : transition){
            for(int i = 0; i < 2; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
            }
        }
    }
    else {
        currentQMatrix = Matrix<double>(61, 61, 0.0);
        for(auto coord : valid){
            currentQMatrix(coord.first, coord.second) = currentStationary[coord.second];
            currentQMatrix(coord.second, coord.first) = currentStationary[coord.first];
        }
        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) *= currentK;
            currentQMatrix(coord.second, coord.first) *= currentK;
        }
    }
}

void RJDPPMatrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;
    oldR = currentR;
    oldRPrior = currentRPrior;
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
        rAcceptCount += 1;
    }

    moveChoice = -1;
}

void RJDPPMatrix::reject() {
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentR = oldR;
    currentRPrior = oldRPrior;
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double RJDPPMatrix::lnPrior() {
    int stateSpace = currentQMatrix.dim1();
    double prior = currentKPrior + currentStationaryPrior;
    if(stateSpace != 61)
        prior += currentRPrior;
    return prior;
}

double RJDPPMatrix::updateK() {
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

    return hastings;
}

double RJDPPMatrix::updateR() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = 2;
    rCount += 1;

    double currentV = currentR;
    double scale = std::exp(rDelta * (rng.uniformRv() - 0.5));
    double newV = currentR * scale;

    currentR = newV;
    hastings = std::log(scale);

    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);

    return hastings;
}

double RJDPPMatrix::updateStationary() {
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

    currentStationaryPrior = Probability::Dirichlet::lnPdf(stationaryPriorAlpha, currentStationary);

    return hastings;
}

Matrix<double> RJDPPMatrix::Q(double omega1, double omega2, double omega3) {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    int stateSpace = currentQMatrix.dim1();

    if(stateSpace == 183){
        for(auto coord : nonsynonymous){
            returnMatrix(coord.first, coord.second) *= omega1;
            returnMatrix(coord.second, coord.first) *= omega1;

            returnMatrix(coord.first + 61, coord.second + 61) *= omega2;
            returnMatrix(coord.second + 61, coord.first + 61) *= omega2;

            returnMatrix(coord.first + 122, coord.second + 122) *= omega3;
            returnMatrix(coord.second + 122, coord.first + 122) *= omega3;
        }
    }
    else if(stateSpace == 122){
        for(auto coord : nonsynonymous){
            returnMatrix(coord.first, coord.second) *= omega1;
            returnMatrix(coord.second, coord.first) *= omega1;

            returnMatrix(coord.first + 61, coord.second + 61) *= omega2;
            returnMatrix(coord.second + 61, coord.first + 61) *= omega2;
        }
    }
    else {
        for(auto coord : nonsynonymous){
            returnMatrix(coord.first, coord.second) *= omega1;
            returnMatrix(coord.second, coord.first) *= omega1;
        }
    }

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
    
    if(stateSpace == 183){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = currentR;
            returnMatrix(i + 61, i) = currentR;
            returnMatrix(i, i + 122) = currentR;
            returnMatrix(i + 122, i) = currentR;
            returnMatrix(i + 61, i + 122) = currentR;
            returnMatrix(i + 122, i + 61) = currentR;

            returnMatrix(i, i) -= 2.0*currentR;
            returnMatrix(i + 61, i + 61) -= 2.0*currentR;
            returnMatrix(i + 122, i + 122) -= 2.0*currentR;
        }
    }
    else if(stateSpace == 122){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = currentR;
            returnMatrix(i + 61, i) = currentR;

            returnMatrix(i, i) -= currentR;
            returnMatrix(i + 61, i + 61) -= currentR;
        } 
    }

    return returnMatrix;
}

std::vector<double> RJDPPMatrix::getStationary(int omegaCount){
    std::vector<double> returnStationary;

    if(omegaCount == 3){
        for(int i = 0; i < 3; i++){
            for(double v : currentStationary){
                returnStationary.push_back(v/3.0);
            }
        }
    }
    else if(omegaCount == 2){
        for(int i = 0; i < 2; i++){
            for(double v : currentStationary){
                returnStationary.push_back(v/2.0);
            }
        }
    }
    else if(omegaCount == 1){
        for(double v : currentStationary){
            returnStationary.push_back(v);
        }
    }
    
    return returnStationary;
}

std::tuple<double, double, double> RJDPPMatrix::dNdS(double omega1, double omega2, double omega3){
    Matrix<double> tempMatrix(183, 183, 0.0);

    for(auto coord : valid){
        for(int i = 0; i < 3; i++){
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) = 1.0;
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) = 1.0;
        }
    }
    for(auto coord : nonsynonymous){
        tempMatrix(coord.first, coord.second) *= omega1;
        tempMatrix(coord.second, coord.first) *= omega1;

        tempMatrix(coord.first + 61, coord.second + 61) *= omega2;
        tempMatrix(coord.second + 61, coord.first + 61) *= omega2;

        tempMatrix(coord.first + 122, coord.second + 122) *= omega3;
        tempMatrix(coord.second + 122, coord.first + 122) *= omega3;
    }
    for(auto coord : transition){
        for(int i = 0; i < 3; i++){
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
        }
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

        dS3 += tempMatrix(coord.first + 122, coord.second + 122) * currentStationary[coord.first];
        dS3 += tempMatrix(coord.second + 122, coord.first + 122) * currentStationary[coord.second];
    }

    return std::make_tuple(dN1/dS1, dN2/dS2, dN3/dS3);
}

void RJDPPMatrix::tune(){
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

    if(rCount > 0){
        double rRate = (double)rAcceptCount/(double)rCount;

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
}