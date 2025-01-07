#include "CodonMultiMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

CodonMultiMatrix::CodonMultiMatrix(Settings settings, std::vector<double> pi) : 
                                   currentQMatrix(122, 122, 0.0), oldQMatrix(122, 122, 0.0), 
                                   currentStationary(pi), oldStationary(pi), kLambda(settings.kLambda), rLambda(settings.rLambda),
                                   currentKPrior(0.0), oldKPrior(0.0), updateStationary(settings.updateStationary), 
                                   currentRPrior(0.0), oldRPrior(0.0), moveChoice(-1), kCount(0), stationaryCount(0), 
                                   kAcceptCount(0), rCount(0), rAcceptCount(0), stationaryAcceptCount(0), 
                                   kDelta(0.5), stationaryAlpha(100), stationaryOffset(0.01), rDelta(0.5) {
    
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
                if(isTransition)
                    transition.insert(pair);
            }
        }
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    if(settings.kValue == -1)
        currentK = Probability::Exponential::rv(&rng, kLambda);
    else
        currentK = settings.kValue;
    oldK = currentK;
    currentKPrior = Probability::Exponential::lnPdf(kLambda, currentK);
    oldKPrior = currentKPrior;


    if(settings.rValue == -1)
        currentR = Probability::Exponential::rv(&rng, rLambda);
    else
        currentR = settings.kValue;
    oldR = currentR;
    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);
    oldRPrior = currentRPrior;

    for(auto coord : valid){
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.second]/2;
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.first]/2;
        currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.second]/2;
        currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.first]/2;

        currentQMatrix(coord.first, coord.first + 61) = currentR * currentStationary[coord.first]/2;
        currentQMatrix(coord.first + 61, coord.first) = currentR * currentStationary[coord.first]/2;
        currentQMatrix(coord.second, coord.second + 61) = currentR * currentStationary[coord.second]/2;
        currentQMatrix(coord.second + 61, coord.second) = currentR * currentStationary[coord.second]/2;
    }
    for(auto coord : transition){
        currentQMatrix(coord.first, coord.second) *= currentK;
        currentQMatrix(coord.second, coord.first) *= currentK;
        currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
        currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;  
    }
    
    oldQMatrix = currentQMatrix;

    for(int i = 0; i < 61; i++)
        randomStates.push_back(i);

    dirty();
}

void CodonMultiMatrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;

    oldR = currentR;
    oldRPrior = currentRPrior;

    oldQMatrix = currentQMatrix;

    if(updateStationary){
        oldStationary = currentStationary;
    }

    if(moveChoice != -1){
        if(moveChoice == 0){
            kAcceptCount += 1;
        }
        else if(moveChoice == 1){
            rAcceptCount += 1;
        }
        else if(moveChoice == 2){
            stationaryAcceptCount += 1;
        }
        moveChoice = -1;
    }
}

void CodonMultiMatrix::reject() {
    currentK = oldK;
    currentKPrior = oldKPrior;

    currentRPrior = oldRPrior;
    currentR = oldR;
    
    currentQMatrix = oldQMatrix;

    if(updateStationary){
        currentStationary = oldStationary;
    }

    moveChoice = -1;
}

double CodonMultiMatrix::lnPrior() {
    return currentKPrior + currentRPrior;
}

double CodonMultiMatrix::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    double randomMove = 0.0;
    do {
        randomMove = rng.uniformRv();
    }
    while(updateStationary == false && randomMove >= 0.66);


    if(randomMove < 0.33) { // Scale K
        moveChoice = 0;
        kCount += 1;

        double currentV = currentK;
        double scale = std::exp(kDelta * (rng.uniformRv() - 0.5));
        double newV = currentK * scale;

        currentK = newV;
        hastings = std::log(scale);

        this->dirty();

        currentKPrior = Probability::Exponential::rv(&rng, kLambda);

        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) = currentStationary[coord.second]/2 * currentK;
            currentQMatrix(coord.second, coord.first) = currentStationary[coord.first]/2 * currentK;
            currentQMatrix(coord.first + 61, coord.second + 61) = currentStationary[coord.second]/2 * currentK;
            currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.first]/2 * currentK;  
        }
    }
    else if(randomMove < 0.66) { // Scale R
        moveChoice = 1;
        rCount += 1;

        double currentV = currentR;
        double scale = std::exp(rDelta * (rng.uniformRv() - 0.5));
        double newV = currentV * scale;

        currentR = newV;
        hastings = std::log(scale);

        this->dirty();

        currentRPrior = Probability::Exponential::rv(&rng, rLambda);

        for(auto coord : valid){
            currentQMatrix(coord.first, coord.first + 61) = currentR * currentStationary[coord.first]/2;
            currentQMatrix(coord.first + 61, coord.first) = currentR * currentStationary[coord.first]/2;
            currentQMatrix(coord.second, coord.second + 61) = currentR * currentStationary[coord.second]/2;
            currentQMatrix(coord.second + 61, coord.second) = currentR * currentStationary[coord.second]/2;
        }
    }
    else { // Dirichlet simplex on the stationary
        moveChoice = 2;
        stationaryCount += 1;

        this->dirty();

        std::vector<int> drawSet(randomStates);
        std::vector<int> randomIndices;

        for(int c = 0; c < 30; c++){
            int i = (int)(rng.uniformRv() * drawSet.size());
            randomIndices.push_back(drawSet[i]);
            drawSet.erase(drawSet.begin() + i);
        }

        std::vector<double> x(31, 0.0);
        std::vector<double> alphaForward(31, 0.0);
        std::vector<double> alphaReverse(31, 0.0);
        std::vector<double> z(31, 0.0);


        for(int i = 0; i < 61; i++) {
            auto it = std::find(randomIndices.begin(), randomIndices.end(), i);
            if(it != randomIndices.end()) {
                x[it - randomIndices.begin()] += currentStationary[i];
            }
            else {
                x[30] += currentStationary[i];
            }
        }

        for(int i = 0; i < x.size(); i++) {
            alphaForward[i] = (x[i] + stationaryOffset) * stationaryAlpha;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = (z[i] + stationaryOffset) * stationaryAlpha;
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
        hastings += 30 * log(factor);

        for(auto coord : valid){
            currentQMatrix(coord.first, coord.second) = currentStationary[coord.second]/2;
            currentQMatrix(coord.second, coord.first) = currentStationary[coord.first]/2;
            currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.second]/2;
            currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.first]/2;

            currentQMatrix(coord.first, coord.first + 61) = currentR * currentStationary[coord.first]/2;
            currentQMatrix(coord.first + 61, coord.first) = currentR * currentStationary[coord.first]/2;
            currentQMatrix(coord.second, coord.second + 61) = currentR * currentStationary[coord.second]/2;
            currentQMatrix(coord.second + 61, coord.second) = currentR * currentStationary[coord.second]/2;
        }
        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) *= currentK;
            currentQMatrix(coord.second, coord.first) *= currentK;
            currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
            currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;  
        }
    }

    return hastings;
}

Matrix<double> CodonMultiMatrix::Q(double omega1, double omega2) {
    Matrix<double> returnMatrix = currentQMatrix.copy();

    for(auto coord : nonsynonymous){
        returnMatrix(coord.first, coord.second) *= omega1;
        returnMatrix(coord.second, coord.first) *= omega1; 

        returnMatrix(coord.first + 61, coord.second + 61) *= omega2;
        returnMatrix(coord.second + 61, coord.first + 61) *= omega2; 
    }

    for(int i = 0; i < 122; i++){
        double total = 0.0;
        for(int j = 0; j < 122; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i , i) = total * -1;
    }

    return returnMatrix;
}

std::vector<double> CodonMultiMatrix::getStationary(){
    std::vector<double> retunStationary;

    for(int i = 0; i < 2; i++){
        for(double v : currentStationary){
            retunStationary.push_back(v/2);
        }
    }
    
    return retunStationary;
}

void CodonMultiMatrix::tune(){
    double rRate = (double)rAcceptCount/(double)rCount;
    if ( rRate > 0.44 ) {
        rDelta *= (1.0 + ((rRate-0.44)/0.766));
    }
    else {
        rDelta /= (2.0 - rRate/0.44);
    }
    rAcceptCount = 0;
    rCount = 0;

    double kRate = (double)kAcceptCount/(double)kCount;

    if ( kRate > 0.44 ) {
        kDelta *= (1.0 + ((kRate-0.44)/0.766));
    }
    else {
        kDelta /= (2.0 - kRate/0.44);
    }
    kAcceptCount = 0;
    kCount = 0;

    if(updateStationary){
        double stationaryRate = (double)stationaryAcceptCount/(double)stationaryCount;

        // If we are getting acceptances too often we should decrease alpha and increase the offset
        if ( stationaryRate > 0.44 ) {
            stationaryOffset *= (1.0 + ((stationaryRate-0.44)/0.766));
            stationaryAlpha /= (2.0 - stationaryRate/0.44);
        }
        else {
            stationaryAlpha *= (1.0 + ((stationaryRate-0.44)/0.766));
            stationaryOffset /= (2.0 - stationaryRate/0.44);
        }

        stationaryAcceptCount = 0;
        stationaryCount = 0;
    }
}