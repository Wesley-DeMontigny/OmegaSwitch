#include "CodonMultiMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

CodonMultiMatrix::CodonMultiMatrix(Settings settings) : 
                                   currentQMatrix(122, 122, 0.0), oldQMatrix(122, 122, 0.0), 
                                   currentStationary(61, -1), oldStationary(61, -1), kAlpha(settings.kAlpha), rAlpha(settings.rAlpha),
                                   currentKPrior(0.0), oldKPrior(0.0), currentRPrior(0.0), oldRPrior(0.0), 
                                   moveChoice(-1), kCount(0), stationaryCount(0),
                                   kAcceptCount(0), rCount(0), rAcceptCount(0), stationaryAcceptCount(0),
                                   kDelta(0.5), stationaryAlpha(25), rDelta(0.5) {
    
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
        currentK = Probability::Gamma::rv(&rng, kAlpha, 1);
    else
        currentK = settings.kValue;
    oldK = currentK;
    currentKPrior = Probability::Gamma::lnPdf(kAlpha, 1, currentK);
    oldKPrior = currentKPrior;


    if(settings.rValue == -1)
        currentR = Probability::Gamma::rv(&rng, rAlpha, 1);
    else
        currentR = settings.kValue;
    oldR = currentR;
    currentRPrior = Probability::Gamma::lnPdf(rAlpha, 1, currentR);
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
        currentQMatrix(i, i + 61) = currentR * currentStationary[i]/2;
        currentQMatrix(i + 61, i) = currentR * currentStationary[i]/2;
    }
    
    oldQMatrix = currentQMatrix.copy();

    dirty();
}

void CodonMultiMatrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;

    oldR = currentR;
    oldRPrior = currentRPrior;

    oldQMatrix = currentQMatrix.copy();

    oldStationary = currentStationary;

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
    
    currentQMatrix = oldQMatrix.copy();

    currentStationary = oldStationary;

    moveChoice = -1;
}

double CodonMultiMatrix::lnPrior() {
    return currentKPrior + currentRPrior;
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

    currentKPrior = Probability::Gamma::lnPdf(kAlpha, 1, currentK);

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

    moveChoice = 1;
    rCount += 1;

    double currentV = currentR;
    double scale = std::exp(rDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentR = newV;
    hastings = std::log(scale);

    this->dirty();

    currentRPrior = Probability::Gamma::lnPdf(rAlpha, 1, currentR);

    for(int i = 0; i < 61; i++){
        currentQMatrix(i, i + 61) = currentR * currentStationary[i]/2;
        currentQMatrix(i + 61, i) = currentR * currentStationary[i]/2;
    }

    return hastings;
}

double CodonMultiMatrix::updateStationary() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    //double draw = rng.uniformRv();

    /*
    if(draw < 0.2){
        moveChoice = 2;
        stationaryDirichletCount += 1;

        this->dirty();

        std::vector<double> x(61, 0.0);
        std::vector<double> alphaForward(61, 0.0);
        std::vector<double> alphaReverse(61, 0.0);
        std::vector<double> z(61, 0.0);


        for(int i = 0; i < 61; i++) {
            x[i] += currentStationary[i];
            alphaForward[i] = (x[i] * stationaryDirichletAlpha) + 1;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = (z[i] * stationaryDirichletAlpha) + 1;
        }

        for(int i = 0; i < 61; i++) {
            currentStationary[i] = z[i];
            if(currentStationary[i] < 1E-10) {
                return -1 * INFINITY;
            }
        }

        hastings  = Probability::Dirichlet::lnPdf(alphaReverse, x) - Probability::Dirichlet::lnPdf(alphaForward, z);
    }
    */
    {
        moveChoice = 2;
        stationaryCount += 1;

        this->dirty();

        int i = (int)(rng.uniformRv() * 61);

        double oldVal = currentStationary[i];

        double a = stationaryAlpha + 1.0;
        double b = (stationaryAlpha / oldVal) - a + 2.0;
        double newVal = Probability::Beta::rv(&rng, a, b);

        currentStationary[i] = newVal;

        double scalingFactor = (1.0 - newVal)/(1.0 - oldVal);

        double sum = 0.0;
        double hastings = 0.0;
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
        double newA = stationaryAlpha + 1.0;
        double newB = (stationaryAlpha / newVal) - a + 2.0;
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
        currentQMatrix(i, i + 61) = currentR * currentStationary[i]/2;
        currentQMatrix(i + 61, i) = currentR * currentStationary[i]/2;
    }

    return hastings;
}

Matrix<double> CodonMultiMatrix::Q(double omega1, double omega2, int invariant) {
    Matrix<double> returnMatrix = currentQMatrix.copy();

    for(auto coord : nonsynonymous){
        returnMatrix(coord.first, coord.second) *= omega1;
        returnMatrix(coord.second, coord.first) *= omega1; 

        returnMatrix(coord.first + 61, coord.second + 61) *= omega2;
        returnMatrix(coord.second + 61, coord.first + 61) *= omega2; 
    }

    if(invariant == 1){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = 0;
            returnMatrix(i + 61, i) = 0;
        }
    }

    double scaler = 0.0;
    for(int i = 0; i < 122; i++){
        double total = 0.0;
        for(int j = 0; j < 122; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i , i) = total * -1;
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
    double rRate = (double)rAcceptCount/(double)rCount;
    if ( rRate > 0.33 ) {
        rDelta *= (1.0 + ((rRate-0.33)/0.67));
    }
    else {
        rDelta /= (2.0 - rRate/0.33);
    }
    rAcceptCount = 0;
    rCount = 0;

    double kRate = (double)kAcceptCount/(double)kCount;

    if ( kRate > 0.33 ) {
        kDelta *= (1.0 + ((kRate-0.33)/0.67));
    }
    else {
        kDelta /= (2.0 - kRate/0.33);
    }
    kAcceptCount = 0;
    kCount = 0;

    /*
    double stationaryDirichletRate = (double)stationaryDirichletAcceptCount/(double)stationaryDirichletCount;

    if ( stationaryDirichletRate > 0.33 ) {
        stationaryDirichletAlpha /= (1.0 + ((stationaryDirichletRate-0.33)/0.67));
    }
    else {
        stationaryDirichletAlpha *= (2.0 - stationaryDirichletRate/0.33);
    }

    stationaryDirichletAlpha = std::fmin(200, stationaryDirichletAlpha);

    stationaryDirichletAcceptCount = 0;
    stationaryDirichletCount = 0;
    */

    double stationaryRate = (double)stationaryAcceptCount/(double)stationaryCount;

    if ( stationaryRate > 0.33 ) {
        stationaryAlpha /= (1.0 + ((stationaryRate-0.33)/0.67));
    }
    else {
        stationaryAlpha *= (2.0 - stationaryRate/0.33);
    }

    stationaryAlpha = std::fmin(100, stationaryRate);

    stationaryAcceptCount = 0;
    stationaryCount = 0;
}