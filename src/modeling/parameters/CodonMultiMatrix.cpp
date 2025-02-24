#include "CodonMultiMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

CodonMultiMatrix::CodonMultiMatrix(Settings settings) : 
                                   currentQMatrix(122, 122, 0.0), oldQMatrix(122, 122, 0.0), 
                                   currentStationary(61, -1), oldStationary(61, -1), kLambda(settings.kLambda), rLambda(settings.rLambda),
                                   currentKPrior(0.0), oldKPrior(0.0), currentRPrior(0.0), oldRPrior(0.0), moveChoice(-1), kCount(0),
                                   stationaryDirichletCount(0), stationaryBetaCount(0), rCount(0), kAcceptCount(0), stationaryDirichletAcceptCount(0), 
                                   stationaryBetaAcceptCount(0), rAcceptCount(0), kDelta(0.5), stationaryBetaAlpha(50),
                                   stationaryDirichletAlpha(500), rDelta(0.5) {
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
        currentR = settings.rValue;
    oldR = currentR;
    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);
    oldRPrior = currentRPrior;

    std::vector<double> alpha;
    for(int i = 0; i < 61; i++)
        alpha.push_back(1.0);
    
    Probability::Dirichlet::rv(&rng, alpha, currentStationary);
    oldStationary = currentStationary;

    for(auto coord : valid){
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.second]/2;
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.first]/2;
        currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.second]/2;
        currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.first]/2;
    }
    for(auto coord : transition){
        currentQMatrix(coord.first, coord.second) *= currentK;
        currentQMatrix(coord.second, coord.first) *= currentK;
        currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
        currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;  
    }
    for(int i = 0; i < 61; i++){
        currentQMatrix(i, i + 61) = currentStationary[i]/2 * currentR;
        currentQMatrix(i + 61, i) = currentStationary[i]/2 * currentR;
    }
    
    oldQMatrix = currentQMatrix.copy();

    for(int i = 0; i < 61; i++)
        randomStates.push_back(i);

    dirty();
}

void CodonMultiMatrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;
    oldR = currentR;
    oldRPrior = currentRPrior;

    oldQMatrix = currentQMatrix.copy();

    oldStationary = currentStationary;

    if(moveChoice == 0){
        kAcceptCount += 1;
    }
    else if(moveChoice == 1){
        stationaryDirichletAcceptCount += 1;
    }
    else if(moveChoice == 2){
        rAcceptCount += 1;
    }
    else if(moveChoice == 3){
        stationaryBetaAcceptCount += 1;
    }
    moveChoice = -1;
}

void CodonMultiMatrix::reject() {
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentR = oldR;
    currentRPrior = oldRPrior;

    currentQMatrix = oldQMatrix.copy();

    currentStationary = oldStationary;

    moveChoice = -1;
}

double CodonMultiMatrix::lnPrior() {
    return currentKPrior;
}

double CodonMultiMatrix::updateK() {
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

    for(auto coord : transition){
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.second]/2 * currentK;
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.first]/2 * currentK;
        currentQMatrix(coord.first + 61, coord.second + 61) = currentStationary[coord.second]/2 * currentK;
        currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.first]/2 * currentK;  
    }

    return hastings;
}

double CodonMultiMatrix::updateR() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = 2;
    rCount += 1;

    double currentV = currentR;
    double scale = std::exp(rDelta * (rng.uniformRv() - 0.5));
    double newV = currentR * scale;

    currentR = newV;
    hastings = std::log(scale);

    this->dirty();

    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);

    for(int i = 0; i < 61; i++){
        currentQMatrix(i, i + 61) = currentStationary[i]/2 * currentR;
        currentQMatrix(i + 61, i) = currentStationary[i]/2 * currentR;
    }

    return hastings;
}

double CodonMultiMatrix::updateStationary() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    double choice = rng.uniformRv();

    if(choice > 0.5){
        moveChoice = 1;
        stationaryDirichletCount += 1;
        this->dirty();

        std::vector<int> drawSet(randomStates);
        std::vector<int> randomIndices;

        for(int c = 0; c < 5; c++){
            int i = (int)(rng.uniformRv() * drawSet.size());
            randomIndices.push_back(drawSet[i]);
            drawSet.erase(drawSet.begin() + i);
        }

        std::vector<double> x(6, 0.0);
        std::vector<double> alphaForward(6, 0.0);
        std::vector<double> alphaReverse(6, 0.0);
        std::vector<double> z(6, 0.0);


        for(int i = 0; i < 61; i++) {
            auto it = std::find(randomIndices.begin(), randomIndices.end(), i);
            if(it != randomIndices.end()) {
                x[it - randomIndices.begin()] += currentStationary[i];
            }
            else {
                x[5] += currentStationary[i];
            }
        }

        for(int i = 0; i < x.size(); i++) {
            alphaForward[i] = (x[i] * stationaryDirichletAlpha) + 0.001;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = (z[i] * stationaryDirichletAlpha) + 0.001;
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

            if(currentStationary[i] < 1E-25) {
                return -1 * INFINITY;
            }
        }

        hastings  = Probability::Dirichlet::lnPdf(alphaReverse, x) - Probability::Dirichlet::lnPdf(alphaForward, z);
        hastings += 55 * log(factor);
    }   
    else {
        moveChoice = 3;
        stationaryBetaCount += 1;
        this->dirty();
    
        int i = (int)(rng.uniformRv() * 61);
    
        double oldVal = currentStationary[i];
    
        double a = stationaryBetaAlpha + 1.0;
        double b = (stationaryBetaAlpha / oldVal) - a + 2.0;
        double newVal = Probability::Beta::rv(&rng, a, b);
    
        currentStationary[i] = newVal;
    
        double scalingFactor = (1.0 - newVal)/(1.0 - oldVal);
    
        double sum = 0.0;
        for(int j = 0; j < 61; j++){
            if(j != i)
                currentStationary[j] = currentStationary[j] * scalingFactor;
    
            if(currentStationary[j] < 1e-10)
                return -1.0 * INFINITY;
            
            sum += currentStationary[j];
        }
    
        //Normalize to make sure this doesn't drift from 1.0
        for (int j = 0; j < 61; j++) {
            currentStationary[j] = currentStationary[j]/sum;
        }
    
        // The probability of getting our new value
        double forward = Probability::Beta::lnPdf(a, b, newVal);
        double newA = stationaryBetaAlpha + 1.0;
        double newB = (stationaryBetaAlpha / newVal) - a + 2.0;
        // The probability of getting our old value in the future
        double backward = Probability::Beta::lnPdf(newA, newB, oldVal);
        
        hastings = backward - forward;
        
        hastings += 59 * std::log(scalingFactor) - 60 * std::log(sum);
    }

    for(auto coord : valid){
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.second]/2;
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.first]/2;
        currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.second]/2;
        currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.first]/2;
    }
    for(auto coord : transition){
        currentQMatrix(coord.first, coord.second) *= currentK;
        currentQMatrix(coord.second, coord.first) *= currentK;
        currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
        currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;
    }
    for(int i = 0; i < 61; i++){
        currentQMatrix(i, i + 61) = currentStationary[i]/2 * currentR;
        currentQMatrix(i + 61, i) = currentStationary[i]/2 * currentR;
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

    double scaler = 0.0;
    for(int i = 0; i < 122; i++){
        double total = 0.0;
        for(int j = 0; j < 122; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1;
        scaler += returnMatrix(i, i);
    }
	
    scaler = -1.0 / scaler;
    for (int i = 0; i < 122; i++)
        for (int j = 0; j < 122; j++)
            returnMatrix(i, j) *= scaler;

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
    double kRate = (double)kAcceptCount/(double)kCount;

    if ( kRate > 0.33 ) {
        kDelta *= (1.0 + ((kRate-0.33)/0.67));
    }
    else {
        kDelta /= (2.0 - kRate/0.33);
    }
    kAcceptCount = 0;
    kCount = 0;

    double rRate = (double)rAcceptCount/(double)rCount;

    if ( rRate > 0.33 ) {
        rDelta *= (1.0 + ((rRate-0.33)/0.67));
    }
    else {
        rDelta /= (2.0 - rRate/0.33);
    }
    rAcceptCount = 0;
    rCount = 0;

    double stationaryDRate = (double)stationaryDirichletAcceptCount/(double)stationaryDirichletCount;

    if ( stationaryDRate > 0.33 ) {
        stationaryDirichletAlpha /= (1.0 + ((stationaryDRate-0.33)/0.67));
    }
    else {
        stationaryDirichletAlpha *= (2.0 - stationaryDRate/0.33);
    }

    stationaryDirichletAlpha = std::fmin(2500.0, stationaryDirichletAlpha);

    stationaryDirichletAcceptCount = 0;
    stationaryDirichletCount = 0;

    double stationaryBRate = (double)stationaryBetaAcceptCount/(double)stationaryBetaCount;

    if ( stationaryBRate > 0.33 ) {
        stationaryBetaAlpha /= (1.0 + ((stationaryBRate-0.33)/0.67));
    }
    else {
        stationaryBetaAlpha *= (2.0 - stationaryBRate/0.33);
    }

    stationaryBetaAlpha = std::fmin(300.0, stationaryBetaAlpha);

    stationaryBetaAcceptCount = 0;
    stationaryBetaCount = 0;
}