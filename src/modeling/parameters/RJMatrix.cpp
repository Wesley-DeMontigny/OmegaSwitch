#include "RJMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

RJMatrix::RJMatrix(Settings settings) : 
                                   currentQMatrix(305, 305, 0.0), oldQMatrix(305, 305, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), rLambda(settings.rLambda), omegaLambda(settings.omegaLambda), rDelta(0.5),  
                                   stationaryAlpha(30000), kDelta(0.5), omegaDelta(0.5), stationaryPriorAlpha(61, 2.0) {
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

    currentOmega4 = Probability::Exponential::rv(&rng, omegaLambda);
    oldOmega4 = currentOmega4;
    currentOmega4Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega4);
    oldOmega4Prior = currentOmega4Prior;

    currentOmega5 = Probability::Exponential::rv(&rng, omegaLambda);
    oldOmega5 = currentOmega5;
    currentOmega5Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega5);
    oldOmega5Prior = currentOmega5Prior;

    currentR = Probability::Exponential::rv(&rng, rLambda);
    oldR = currentR;
    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);
    oldRPrior = currentRPrior;

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
    
    if(currentActiveOmegas == 5){
        currentQMatrix = Matrix<double>(305, 305, 0.0);
        for(auto coord : valid){
            for(int i = 0; i < 5; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
            }
        }
        for(auto coord : nonsynonymous){
            currentQMatrix(coord.first, coord.second) *= currentOmega1;
            currentQMatrix(coord.second, coord.first) *= currentOmega1;

            currentQMatrix(coord.first + 61, coord.second + 61) *= currentOmega1 + currentOmega2;
            currentQMatrix(coord.second + 61, coord.first + 61) *= currentOmega1 + currentOmega2;

            currentQMatrix(coord.first + 122, coord.second + 122) *= currentOmega1 + currentOmega2 + currentOmega3;
            currentQMatrix(coord.second + 122, coord.first + 122) *= currentOmega1 + currentOmega2 + currentOmega3;

            currentQMatrix(coord.first + 183, coord.second + 183) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4;
            currentQMatrix(coord.second + 183, coord.first + 183) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4;

            currentQMatrix(coord.first + 244, coord.second + 244) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4 + currentOmega5;
            currentQMatrix(coord.second + 244, coord.first + 244) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4 + currentOmega5;            
        }
        for(auto coord : transition){
            for(int i = 0; i < 5; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
            }
        }
    }
    else if(currentActiveOmegas == 4){
        currentQMatrix = Matrix<double>(244, 244, 0.0);
        for(auto coord : valid){
            for(int i = 0; i < 4; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
            }
        }
        for(auto coord : nonsynonymous){
            currentQMatrix(coord.first, coord.second) *= currentOmega1;
            currentQMatrix(coord.second, coord.first) *= currentOmega1;

            currentQMatrix(coord.first + 61, coord.second + 61) *= currentOmega1 + currentOmega2;
            currentQMatrix(coord.second + 61, coord.first + 61) *= currentOmega1 + currentOmega2;

            currentQMatrix(coord.first + 122, coord.second + 122) *= currentOmega1 + currentOmega2 + currentOmega3;
            currentQMatrix(coord.second + 122, coord.first + 122) *= currentOmega1 + currentOmega2 + currentOmega3;

            currentQMatrix(coord.first + 183, coord.second + 183) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4;
            currentQMatrix(coord.second + 183, coord.first + 183) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4;           
        }
        for(auto coord : transition){
            for(int i = 0; i < 4; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
            }
        }
    }
    else if(currentActiveOmegas == 3){
        currentQMatrix = Matrix<double>(183, 183, 0.0);
        for(auto coord : valid){
            for(int i = 0; i < 3; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
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
        for(auto coord : transition){
            for(int i = 0; i < 3; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
            }
        }
    }
    else if(currentActiveOmegas == 2){
        currentQMatrix = Matrix<double>(122, 122, 0.0);
        for(auto coord : valid){
            for(int i = 0; i < 2; i++){
                currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second];
                currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first];
            }
        }
        for(auto coord : nonsynonymous){
            currentQMatrix(coord.first, coord.second) *= currentOmega1;
            currentQMatrix(coord.second, coord.first) *= currentOmega1;

            currentQMatrix(coord.first + 61, coord.second + 61) *= currentOmega1 + currentOmega2;
            currentQMatrix(coord.second + 61, coord.first + 61) *= currentOmega1 + currentOmega2;
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
        for(auto coord : nonsynonymous){
            currentQMatrix(coord.first, coord.second) *= currentOmega1;
            currentQMatrix(coord.second, coord.first) *= currentOmega1;
        }
        for(auto coord : transition){
            currentQMatrix(coord.first, coord.second) *= currentK;
            currentQMatrix(coord.second, coord.first) *= currentK;
        }
    }
}

void RJMatrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;
    oldOmega1 = currentOmega1;
    oldOmega1Prior = currentOmega1Prior;
    oldOmega2 = currentOmega2;
    oldOmega2Prior = currentOmega2Prior;
    oldOmega3 = currentOmega3;
    oldOmega3Prior = currentOmega3Prior;
    oldOmega4 = currentOmega4;
    oldOmega4Prior = currentOmega4Prior;
    oldOmega5 = currentOmega5;
    oldOmega5Prior = currentOmega5Prior;
    oldR = currentR;
    oldRPrior = currentRPrior;
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
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentOmega1 = oldOmega1;
    currentOmega1Prior = oldOmega1Prior;
    currentOmega2 = oldOmega2;
    currentOmega2Prior = oldOmega2Prior;
    currentOmega3 = oldOmega3;
    currentOmega3Prior = oldOmega3Prior;
    currentOmega4 = oldOmega4;
    currentOmega4Prior = oldOmega4Prior;
    currentOmega5 = oldOmega5;
    currentOmega5Prior = oldOmega5Prior;
    currentR = oldR;
    currentRPrior = oldRPrior;
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;
    currentActiveOmegas = oldActiveOmegas;

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double RJMatrix::lnPrior() {
    double prior = currentKPrior + currentStationaryPrior;

    if(currentActiveOmegas == 5){
        prior += currentOmega1Prior + currentOmega2Prior + currentOmega3Prior + currentOmega4Prior + currentOmega5Prior + currentRPrior;
    }
    else if(currentActiveOmegas == 4){
        prior += currentOmega1Prior + currentOmega2Prior + currentOmega3Prior + currentOmega4Prior + currentRPrior;
    }
    else if(currentActiveOmegas == 3){
        prior += currentOmega1Prior + currentOmega2Prior + currentOmega3Prior + currentRPrior;
    }
    else if(currentActiveOmegas == 2){
        prior += currentOmega1Prior + currentOmega2Prior + currentRPrior;
    }
    else if(currentActiveOmegas == 1){
        prior += currentOmega1Prior;
    }

    return prior;
}

double RJMatrix::updateK() {
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

double RJMatrix::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = 2;
    omegaCount += 1;

    double randomOmega = rng.uniformRv();

    if(currentActiveOmegas == 5){
        if(randomOmega < 0.2){
        double currentV = currentOmega1;
        double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV = currentOmega1 * scale;

        currentOmega1 = newV;
        hastings = std::log(scale);

        currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
        }
        else if(randomOmega < 0.4){
            double currentV = currentOmega2;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega2 * scale;

            currentOmega2 = newV;
            hastings = std::log(scale);

            currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);
        }
        else if(randomOmega < 0.6){
            double currentV = currentOmega3;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega3 * scale;

            currentOmega3 = newV;
            hastings = std::log(scale);

            currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);
        }
        else if(randomOmega < 0.8){
            double currentV = currentOmega4;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega4 * scale;

            currentOmega4 = newV;
            hastings = std::log(scale);

            currentOmega4Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega4);
        }
        else{
            double currentV = currentOmega5;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega5 * scale;

            currentOmega5 = newV;
            hastings = std::log(scale);

            currentOmega5Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega5);
        }
    }
    else if(currentActiveOmegas == 4){
        if(randomOmega < 0.25){
            double currentV = currentOmega1;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega1 * scale;

            currentOmega1 = newV;
            hastings = std::log(scale);

            currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
        }
        else if(randomOmega < 0.50){
            double currentV = currentOmega2;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega2 * scale;

            currentOmega2 = newV;
            hastings = std::log(scale);

            currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);
        }
        else if(randomOmega < 0.75){
            double currentV = currentOmega3;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega3 * scale;

            currentOmega3 = newV;
            hastings = std::log(scale);

            currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);
        }
        else{
            double currentV = currentOmega4;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega4 * scale;

            currentOmega4 = newV;
            hastings = std::log(scale);

            currentOmega4Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega4);
        }
    }
    else if(currentActiveOmegas == 3){
        if(randomOmega < 0.33){
            double currentV = currentOmega1;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega1 * scale;

            currentOmega1 = newV;
            hastings = std::log(scale);

            currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
        }
        else if(randomOmega < 0.66){
            double currentV = currentOmega2;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega2 * scale;

            currentOmega2 = newV;
            hastings = std::log(scale);

            currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);
        }
        else{
            double currentV = currentOmega3;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega3 * scale;

            currentOmega3 = newV;
            hastings = std::log(scale);

            currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);
        }
    }
    else if(currentActiveOmegas == 2){
        if(randomOmega < 0.50){
            double currentV = currentOmega1;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega1 * scale;

            currentOmega1 = newV;
            hastings = std::log(scale);

            currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
        }
        else{
            double currentV = currentOmega2;
            double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV = currentOmega2 * scale;

            currentOmega2 = newV;
            hastings = std::log(scale);

            currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);
        }
    }
    else {
        double currentV = currentOmega1;
        double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV = currentOmega1 * scale;

        currentOmega1 = newV;
        hastings = std::log(scale);

        currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
    }

    rebuildQMatrix();

    return hastings;
}

double RJMatrix::updateR() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = 3;
    rCount += 1;

    double currentV = currentR;
    double scale = std::exp(rDelta * (rng.uniformRv() - 0.5));
    double newV = currentR * scale;

    currentR = newV;
    hastings = std::log(scale);

    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);

    return hastings;
}

double RJMatrix::updateActiveOmegas() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double alpha = 15.0;
    double hastings = 0.0;
    this->dirty();

    if (currentActiveOmegas == 1) { // 1 2 split
        hastings = std::log(2.0);
        currentActiveOmegas = 2;

        double tempO = currentOmega1;
        double u = Probability::Beta::rv(&rng, alpha, alpha);
        currentOmega1 = tempO * u;
        currentOmega2 = tempO * (1.0 - u);
        currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);
        currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);

        currentR = Probability::Exponential::rv(&rng, rLambda);
        currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);

        hastings += std::log(tempO);
        hastings += Probability::Beta::lnPdf(alpha, alpha, u) + currentRPrior;
    }
    else if (currentActiveOmegas == 2) { // 2 1 merge or 2 3 split
        hastings = std::log(0.5);

        if (rng.uniformRv() > 0.5) { // 2 3 split
            currentActiveOmegas = 3;

            double tempO = currentOmega2;
            double u = Probability::Beta::rv(&rng, alpha, alpha);
            currentOmega2 = tempO * u;
            currentOmega3 = tempO * (1.0 - u);
            currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);
            currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);

            hastings += std::log(tempO);
            hastings += Probability::Beta::lnPdf(alpha, alpha, u);
        } else { // 2 1 merge
            currentActiveOmegas = 1;
            double o1 = currentOmega1;
            double o2 = currentOmega2;
            double total = o1 + o2;
            double u = o1 / total;

            currentOmega1 = total;
            currentOmega1Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega1);

            hastings -= std::log(total);
            hastings -= Probability::Beta::lnPdf(alpha, alpha, u) + Probability::Exponential::lnPdf(rLambda, currentR);
        }
    }
    else if (currentActiveOmegas == 3) { // 3 2 merge or 3 4 split
        hastings = std::log(0.5);

        if (rng.uniformRv() > 0.5) { // 3 4 split
            currentActiveOmegas = 4;

            double tempO = currentOmega3;
            double u = Probability::Beta::rv(&rng, alpha, alpha);
            currentOmega3 = tempO * u;
            currentOmega4 = tempO * (1.0 - u);
            currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);
            currentOmega4Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega4);

            hastings += std::log(tempO);
            hastings += Probability::Beta::lnPdf(alpha, alpha, u);
        } else { // 3 2 merge
            currentActiveOmegas = 2;
            double o2 = currentOmega2;
            double o3 = currentOmega3;
            double total = o2 + o3;
            double u = o2 / total;

            currentOmega2 = total;
            currentOmega2Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega2);

            hastings -= std::log(total);
            hastings -= Probability::Beta::lnPdf(alpha, alpha, u);
        }
    }
    else if (currentActiveOmegas == 4) { // 4 3 merge or 4 5 split
        hastings = std::log(0.5);

        if (rng.uniformRv() > 0.5) { // 4 5 split
            currentActiveOmegas = 5;

            double tempO = currentOmega4;
            double u = Probability::Beta::rv(&rng, alpha, alpha);
            currentOmega4 = tempO * u;
            currentOmega5 = tempO * (1.0 - u);
            currentOmega4Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega4);
            currentOmega5Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega5);

            hastings += std::log(tempO);
            hastings += Probability::Beta::lnPdf(alpha, alpha, u);
        } else { // 4 3 merge
            currentActiveOmegas = 3;
            double o3 = currentOmega3;
            double o4 = currentOmega4;
            double total = o3 + o4;
            double u = o3 / total;

            currentOmega3 = total;
            currentOmega3Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega3);

            hastings -= std::log(total);
            hastings -= Probability::Beta::lnPdf(alpha, alpha, u);
        }
    }
    else { // 5 4 merge
        currentActiveOmegas = 4;
        hastings = std::log(0.5);

        double o4 = currentOmega4;
        double o5 = currentOmega5;
        double total = o4 + o5;
        double u = o4 / total;

        currentOmega4 = total;
        currentOmega4Prior = Probability::Exponential::lnPdf(omegaLambda, currentOmega4);

        hastings -= std::log(total);
        hastings -= Probability::Beta::lnPdf(alpha, alpha, u);
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
                    returnMatrix(i + offsets[a], i + offsets[b]) = currentR;
                    returnMatrix(i + offsets[b], i + offsets[a]) = currentR;
                    returnMatrix(i + offsets[a], i + offsets[a]) -= currentR;
                    returnMatrix(i + offsets[b], i + offsets[b]) -= currentR;
                }
            }
        }
    }
    else if(currentActiveOmegas == 4){
        for(int i = 0; i < 61; i++){
            int offsets[4] = {0, 61, 122, 183};
            for(int a = 0; a < 4; a++){
                for(int b = a + 1; b < 4; b++){
                    returnMatrix(i + offsets[a], i + offsets[b]) = currentR;
                    returnMatrix(i + offsets[b], i + offsets[a]) = currentR;
                    returnMatrix(i + offsets[a], i + offsets[a]) -= currentR;
                    returnMatrix(i + offsets[b], i + offsets[b]) -= currentR;
                }
            }
        }
    }
    else if(currentActiveOmegas == 3){
        for(int i = 0; i < 61; i++){
            int offsets[3] = {0, 61, 122};
            for(int a = 0; a < 3; a++){
                for(int b = a + 1; b < 3; b++){
                    returnMatrix(i + offsets[a], i + offsets[b]) = currentR;
                    returnMatrix(i + offsets[b], i + offsets[a]) = currentR;
                    returnMatrix(i + offsets[a], i + offsets[a]) -= currentR;
                    returnMatrix(i + offsets[b], i + offsets[b]) -= currentR;
                }
            }
        }
    }
    else if(currentActiveOmegas == 2){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = currentR;
            returnMatrix(i + 61, i) = currentR;
            returnMatrix(i, i) -= currentR;
            returnMatrix(i + 61, i + 61) -= currentR;
        }
    }


    return returnMatrix;
}

std::vector<double> RJMatrix::getStationary(){
    std::vector<double> returnStationary;

    for(int i = 0; i < currentActiveOmegas; i++){
        for(double v : currentStationary){
            returnStationary.push_back(v/(double)(currentActiveOmegas));
        }
    }
    
    return returnStationary;
}

std::tuple<double, double, double, double, double> RJMatrix::dNdS(){
    Matrix<double> tempMatrix(305, 305, 0.0);

    for(auto coord : valid){
        for(int i = 0; i < 5; i++){
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) = 1.0;
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) = 1.0;
        }
    }
    for(auto coord : nonsynonymous){
        tempMatrix(coord.first, coord.second) *= currentOmega1;
        tempMatrix(coord.second, coord.first) *= currentOmega1;

        tempMatrix(coord.first + 61, coord.second + 61) *= currentOmega1 + currentOmega2;
        tempMatrix(coord.second + 61, coord.first + 61) *= currentOmega1 + currentOmega2;

        tempMatrix(coord.first + 122, coord.second + 122) *= currentOmega1 + currentOmega2 + currentOmega3;
        tempMatrix(coord.second + 122, coord.first + 122) *= currentOmega1 + currentOmega2 + currentOmega3;

        tempMatrix(coord.first + 183, coord.second + 183) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4;
        tempMatrix(coord.second + 183, coord.first + 183) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4;

        tempMatrix(coord.first + 244, coord.second + 244) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4 + currentOmega5;
        tempMatrix(coord.second + 244, coord.first + 244) *= currentOmega1 + currentOmega2 + currentOmega3 + currentOmega4 + currentOmega5;
    }
    for(auto coord : transition){
        for(int i = 0; i < 5; i++){
            tempMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
            tempMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
        }
    }

    double dN1 = 0.0;
    double dN2 = 0.0;
    double dN3 = 0.0;
    double dN4 = 0.0;
    double dN5 = 0.0;
    for(auto coord : nonsynonymous){
        dN1 += tempMatrix(coord.first, coord.second) * currentStationary[coord.first];
        dN1 += tempMatrix(coord.second, coord.first) * currentStationary[coord.second];

        dN2 += tempMatrix(coord.first + 61, coord.second + 61) * currentStationary[coord.first];
        dN2 += tempMatrix(coord.second + 61, coord.first + 61) * currentStationary[coord.second];

        dN3 += tempMatrix(coord.first + 122, coord.second + 122) * currentStationary[coord.first];
        dN3 += tempMatrix(coord.second + 122, coord.first + 122) * currentStationary[coord.second];

        dN4 += tempMatrix(coord.first + 183, coord.second + 183) * currentStationary[coord.first];
        dN4 += tempMatrix(coord.second + 183, coord.first + 183) * currentStationary[coord.second];

        dN5 += tempMatrix(coord.first + 244, coord.second + 244) * currentStationary[coord.first];
        dN5 += tempMatrix(coord.second + 244, coord.first + 244) * currentStationary[coord.second];
    }

    double dS1 = 0.0;
    double dS2 = 0.0;
    double dS3 = 0.0;
    double dS4 = 0.0;
    double dS5 = 0.0;
    for(auto coord : synonymous){
        dS1 += tempMatrix(coord.first, coord.second) * currentStationary[coord.first];
        dS1 += tempMatrix(coord.second, coord.first) * currentStationary[coord.second];

        dS2 += tempMatrix(coord.first + 61, coord.second + 61) * currentStationary[coord.first];
        dS2 += tempMatrix(coord.second + 61, coord.first + 61) * currentStationary[coord.second];

        dS3 += tempMatrix(coord.first + 122, coord.second + 122) * currentStationary[coord.first];
        dS3 += tempMatrix(coord.second + 122, coord.first + 122) * currentStationary[coord.second];
    
        dS4 += tempMatrix(coord.first + 183, coord.second + 183) * currentStationary[coord.first];
        dS4 += tempMatrix(coord.second + 183, coord.first + 183) * currentStationary[coord.second];

        dS5 += tempMatrix(coord.first + 244, coord.second + 244) * currentStationary[coord.first];
        dS5 += tempMatrix(coord.second + 244, coord.first + 244) * currentStationary[coord.second];
    }

    return std::make_tuple(dN1/dS1, dN2/dS2, dN3/dS3, dN4/dS4, dN5/dS5);
}

void RJMatrix::tune(){
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