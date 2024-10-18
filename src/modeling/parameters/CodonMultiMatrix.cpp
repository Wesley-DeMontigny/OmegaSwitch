#include "CodonMultiMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"

CodonMultiMatrix::CodonMultiMatrix(double rL, std::vector<double> pi, bool updatePi) : 
                                   currentQMatrix(122, 122, 0.0), oldQMatrix(122, 122, 0.0), 
                                   currentStationary(pi), oldStationary(pi), currentRPrior(0.0), oldRPrior(0.0),
                                   rLambda(rL), currentKPrior(0.0), oldKPrior(0.0), currentStationaryPrior(0.0), 
                                   oldStationaryPrior(0.0), updateStationary(updatePi) {
    
    std::vector<int> aaMap = {8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 19, 19, 15, 15, 15, 15, 1, 18, 1, 9, 4, 9, 4};  
    std::vector<char*> codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

    // Because of the complicated nature of this matrix, we need to classify each of the positions in the matrix;
    for(int i = 0; i < 61; i++){
        for(int j = i + 1; j < 61; j++){
            int mismatch = 0;
            for(int k = 0; k < 3; k++){
                if(codons[i][k] != codons[j][k]){
                    mismatch++;
                    if(mismatch > 1){
                        break;
                    }
                    if((codons[i][k] == 'A' && codons[j][k] == 'G') || (codons[i][k] == 'G' && codons[j][k] == 'A') || 
                       (codons[i][k] == 'T' && codons[j][k] == 'C') || (codons[i][k] == 'C' && codons[j][k] == 'T'))
                        transition.insert(std::make_pair(i, j));
                }
            }
            if(mismatch == 1){
                valid.insert(std::make_pair(i, j));
                if(aaMap[i] != aaMap[j])
                    nonsynonymous.insert(std::make_pair(i, j));
            }
        }
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    currentR = Probability::Exponential::rv(&rng, rLambda);
    oldR = currentR;
    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);
    oldRPrior = currentRPrior;

    currentK = (Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1));
    oldK = currentK;
    currentKPrior = -2 * std::log(1.0 + currentK);
    oldKPrior = currentKPrior;

    if(updateStationary){
        for(int i = 0; i < 122; i++){
            flatDirichlet.push_back(1.0);
        }
        currentStationaryPrior = Probability::Dirichlet::lnPdf(flatDirichlet, currentStationary);
        oldStationaryPrior = currentStationaryPrior;
    }

    for(auto coord : valid){
        currentQMatrix(coord.first, coord.second) = currentStationary[coord.first];
        currentQMatrix(coord.second, coord.first) = currentStationary[coord.second];
        currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.first];
        currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.second];

        currentQMatrix(coord.first, coord.first + 61) = currentR * currentStationary[coord.first];
        currentQMatrix(coord.first + 61, coord.first) = currentR * currentStationary[coord.first];
        currentQMatrix(coord.second, coord.second + 61) = currentR * currentStationary[coord.second];
        currentQMatrix(coord.second + 61, coord.second) = currentR * currentStationary[coord.second];
    }
    for(auto coord : transition){
        currentQMatrix(coord.first, coord.second) *= currentK;
        currentQMatrix(coord.second, coord.first) *= currentK;

        currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
        currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;  
    }
    
    oldQMatrix = currentQMatrix;
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
        oldStationaryPrior = currentStationaryPrior;
    }
}

void CodonMultiMatrix::reject() {
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentR = oldR;
    currentRPrior = oldRPrior;

    currentQMatrix = oldQMatrix;

    if(updateStationary){
        currentStationary = oldStationary;
        currentStationaryPrior = oldStationaryPrior;
    }
}

double CodonMultiMatrix::lnPrior() {
    if(!updateStationary) {
        return currentKPrior + currentRPrior;
    }
    else {
        return currentKPrior + currentRPrior + currentStationaryPrior;
    }
}

double CodonMultiMatrix::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    double randomMove = 0.0;
    do {
        randomMove = rng.uniformRv();
    }
    while(updateStationary == false && randomMove >= 0.66);

    if(randomMove < 0.33) { // Scale R
        double delta = std::log(4);

        double scale = std::exp(delta * (rng.uniformRv() - 0.5));
        currentR *= scale;

        this->dirty();

        currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);

        hastings = std::log(scale);

        for(auto coord : valid){
            currentQMatrix(coord.first, coord.first + 61) = currentR * currentStationary[coord.first];
            currentQMatrix(coord.first + 61, coord.first) = currentR * currentStationary[coord.first];
            currentQMatrix(coord.second, coord.second + 61) = currentR * currentStationary[coord.second];
            currentQMatrix(coord.second + 61, coord.second) = currentR * currentStationary[coord.second];
        }

    }
    else if(randomMove < 0.66) { // Scale K
        double delta = std::log(4);

        double scale = std::exp(delta * (rng.uniformRv() - 0.5));
        currentK *= scale;

        this->dirty();

        currentKPrior = -2 * std::log(1.0 + currentK);

        hastings = std::log(scale);

        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) = currentStationary[coord.first];
            currentQMatrix(coord.second, coord.first) = currentStationary[coord.second];
            currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.first];
            currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.second];

            currentQMatrix(coord.first, coord.second) *= currentK;
            currentQMatrix(coord.second, coord.first) *= currentK;

            currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
            currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;  
        }

        currentKPrior = -2 * std::log(1.0 + currentK);

    }
    else { // Beta simplex on the stationary
        int paramNum = currentStationary.size();
        int i = (int)(rng.uniformRv() * paramNum);
        double alpha = 1.0;

        double pVal = currentStationary[i];

        double a = alpha + 1.0;
        double b = (alpha / pVal) - a + 2.0;
        double newVal = Probability::Beta::rv(&rng, a, b);

        currentStationary[i] = newVal;

        double scalingFactor = (1.0 - newVal)/(1.0 - pVal);

        double sum = 0.0;
        for(int j = 0; j < paramNum; j++){
            if(j != i)
                currentStationary[j] *= scalingFactor;

            if(currentStationary[j] <= 1E-100)
                return -1.0 * INFINITY;
            
            sum += currentStationary[j];
        }

        //Normalize to make sure this doesn't drift from 1.0;
        for (int j = 0; j < paramNum; j++) {
            currentStationary[j] /= sum;
        }

        // The probability of getting our new value
        double forward = Probability::Beta::lnPdf(a, b, newVal);
        double newA = alpha + 1.0;
        double newB = (alpha / newVal) - a + 2.0;
        // The probability of getting our old value in the future
        double backward = Probability::Beta::lnPdf(newA, newB, pVal);
        
        hastings = backward - forward;
        
        hastings += (paramNum - 2) * std::log(scalingFactor) - (paramNum - 1) * std::log(sum);

        this->dirty();

        for(auto coord : valid){
            currentQMatrix(coord.first, coord.second) = currentStationary[coord.first];
            currentQMatrix(coord.second, coord.first) = currentStationary[coord.second];
            currentQMatrix(coord.first + 61, coord.second +  61) = currentStationary[coord.first];
            currentQMatrix(coord.second + 61, coord.first + 61) = currentStationary[coord.second];

            currentQMatrix(coord.first, coord.first + 61) = currentR * currentStationary[coord.first];
            currentQMatrix(coord.first + 61, coord.first) = currentR * currentStationary[coord.first];
            currentQMatrix(coord.second, coord.second + 61) = currentR * currentStationary[coord.second];
            currentQMatrix(coord.second + 61, coord.second) = currentR * currentStationary[coord.second];
        }
        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) *= currentK;
            currentQMatrix(coord.second, coord.first) *= currentK;

            currentQMatrix(coord.first + 61, coord.second + 61) *= currentK;
            currentQMatrix(coord.second + 61, coord.first + 61) *= currentK;  
        }

        currentStationaryPrior = Probability::Dirichlet::lnPdf(flatDirichlet, currentStationary);
    }

    return hastings;
}

Matrix<double> CodonMultiMatrix::Q(double omega1, double omega2) {
    Matrix<double> returnMatrix = currentQMatrix;

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
