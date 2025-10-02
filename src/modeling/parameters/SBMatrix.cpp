#include "SBMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

SBMatrix::SBMatrix(Settings settings) : 
                                   currentQMatrix(61 * settings.truncation, 61 * settings.truncation, 0.0), oldQMatrix(61 * settings.truncation, 61 * settings.truncation, 0.0), 
                                   currentStationary(61, -1), oldStationary(61, -1), kLambda(settings.kLambda), rLambda(settings.rLambda),
                                   omegaLambda(settings.omegaLambda), rDelta(0.5),  stationaryAlpha(30000), kDelta(0.5), omegaDelta(0.5),
                                   stationaryPriorAlpha(61, 2.0), proportionAlpha(10), truncation(settings.truncation) {
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

    for(int i = 0; i < truncation; i++){
        double newOmega = Probability::Exponential::rv(&rng, omegaLambda);
        currentOmegas.push_back(newOmega);
        currentOmegasPrior += Probability::Exponential::lnPdf(omegaLambda, newOmega);
    }
    oldOmegasPrior = currentOmegasPrior;
    oldOmegas = currentOmegas;

    currentR = Probability::Exponential::rv(&rng, rLambda);
    oldR = currentR;
    currentRPrior = Probability::Exponential::lnPdf(rLambda, currentR);
    oldRPrior = currentRPrior;
    
    Probability::Dirichlet::rv(&rng, stationaryPriorAlpha, currentStationary);
    oldStationary = currentStationary;
    currentStationaryPrior = Probability::Dirichlet::lnPdf(stationaryPriorAlpha, currentStationary);
    oldStationaryPrior = currentStationaryPrior;

    for(int i = 0; i < truncation -1; i++){
        double newProp = Probability::Beta::rv(&rng, 1.0, 1.0);
        //currentProportionsPrior += Probability::Beta::lnPdf(1.0, 1.0, newProp);
        currentProportionParams.push_back(newProp);
        oldProportionParams.push_back(newProp);
    }

    rebuildQMatrix();
    
    oldQMatrix = currentQMatrix.copy();

    for(int i = 0; i < 61; i++)
        randomStates.push_back(i);

    dirty();
}

void SBMatrix::rebuildQMatrix(){
    proportions.clear();
    proportions.push_back(currentProportionParams[0]);
    double sum = proportions[0];
    for(int i = 1; i < truncation-1; i++){
        double newProp = (1-sum) * currentProportionParams[i];
        proportions.push_back(newProp);
        sum += newProp;
    }
    proportions.push_back(1.0-sum);

    for(auto coord : valid){
        for(int i = 0; i < truncation; i++){
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) = currentStationary[coord.second] * proportions[i];
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) = currentStationary[coord.first] * proportions[i];
        }
    }
    for(auto coord : transition){
        for(int i = 0; i < truncation; i++){
            currentQMatrix(coord.first + (i*61), coord.second + (i*61)) *= currentK;
            currentQMatrix(coord.second + (i*61), coord.first + (i*61)) *= currentK;
        }
    }
    for(auto coord : nonsynonymous){
        for(int i = 0; i < truncation; i++){
            currentQMatrix(coord.first, coord.second) *= currentOmegas[i];
            currentQMatrix(coord.second, coord.first) *= currentOmegas[i];
        }
    }
}

void SBMatrix::accept() {
    oldK = currentK;
    oldKPrior = currentKPrior;
    oldOmegas = currentOmegas;
    oldOmegasPrior = currentOmegasPrior;
    oldR = currentR;
    oldRPrior = currentRPrior;
    oldStationary = currentStationary;
    oldStationaryPrior = currentStationaryPrior;
    oldProportionParams = currentProportionParams;

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
    else if(moveChoice == 4){
        proportionAcceptCount += 1;
    }

    moveChoice = -1;
}

void SBMatrix::reject() {
    currentK = oldK;
    currentKPrior = oldKPrior;
    currentOmegas = oldOmegas;
    currentOmegasPrior = oldOmegasPrior;
    currentR = oldR;
    currentRPrior = oldRPrior;
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;
    currentProportionParams = oldProportionParams;

    proportions.clear();
    proportions.push_back(currentProportionParams[0]);
    double total = (1-proportions[0]);
    for(int i = 1; i < truncation-1; i++){
        double sum = 0.0;
        double newProp = total * proportions[i];
        total *= (1-newProp);
        proportions.push_back(newProp);
    }
    proportions.push_back(1.0-total);

    currentQMatrix = oldQMatrix.copy();

    moveChoice = -1;
}

double SBMatrix::lnPrior() {
    return currentKPrior + currentStationaryPrior + currentOmegasPrior + currentRPrior;
}

double SBMatrix::updateK() {
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

double SBMatrix::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = 2;
    omegaCount += 1;

    double randomOmega = (int)(rng.uniformRv() * truncation);

    double currentV = currentOmegas[randomOmega];
    double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentOmegas[randomOmega] = newV;
    hastings = std::log(scale);


    currentOmegasPrior = 0.0;
    for(int i = 0; i < truncation; i++){
        currentOmegasPrior += Probability::Exponential::lnPdf(omegaLambda, currentOmegas[i]);
    }
    

    rebuildQMatrix();

    return hastings;
}

double SBMatrix::updateR() {
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

double SBMatrix::updateStationary() {
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

double SBMatrix::updateProportions(){
        RandomVariable& rng = RandomVariable::randomVariableInstance();

    moveChoice = 4;
    proportionCount += 1;

    this->dirty();

    int i = (int)(rng.uniformRv() * (truncation-1));

    double oldVal = currentProportionParams[i];

    double a = proportionAlpha+ 1.0;
    double b = (proportionAlpha / oldVal) - a + 2.0;
    double newVal = Probability::Beta::rv(&rng, a, b);

    currentProportionParams[i] = newVal;

    double hastings = 0.0;
    if(currentProportionParams[i] < 1e-5)
        return -1.0 * INFINITY;


    double forward = Probability::Beta::lnPdf(a, b, newVal);
    double newA = proportionAlpha + 1.0;
    double newB = (proportionAlpha / newVal) - a + 2.0;

    double backward = Probability::Beta::lnPdf(newA, newB, oldVal);
    
    hastings = backward - forward;

    rebuildQMatrix();

    return hastings;
}

Matrix<double> SBMatrix::Q() {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    double scaler= 0.0;
    for(int i = 0; i < currentQMatrix.dim1(); i++){
        double total = 0.0;
        for(int j = 0; j < currentQMatrix.dim1(); j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1;
        scaler += returnMatrix(i, i) * currentStationary[i % 61];
    }

    scaler = -1.0 / scaler;
    for (int i = 0; i < currentQMatrix.dim1(); i++)
        for (int j = 0; j < currentQMatrix.dim1(); j++)
            returnMatrix(i, j) *= scaler;
    

    for (int i = 0; i < 61; ++i) {
        for (int c1 = 0; c1 < truncation; ++c1) {
            int index1 = i + c1 * 61;

            for (int c2 = 0; c2 < truncation; ++c2) {
                if (c1 == c2) continue;

                int index2 = i + c2 * 61;
                returnMatrix(index1, index2) = currentR * proportions[c2];
            }

            returnMatrix(index1, index1) -= currentR * (1.0 - proportions[c1]);
        }
    }


    return returnMatrix;
}

std::vector<double> SBMatrix::getStationary(){
    std::vector<double> returnStationary;

    for(int i = 0; i < truncation; i++){
        for(double v : currentStationary){
            returnStationary.push_back(v * proportions[i]);
        }
    }
    
    return returnStationary;
}

std::vector<double> SBMatrix::dNdS(){
        Matrix<double> tempMatrix(currentQMatrix.copy());

    for(auto coord : valid){
        for(int i = 0; i < truncation; i++){
            tempMatrix(coord.first + (61*i), coord.second +  (61*i)) /= currentStationary[coord.second] * proportions[i];
            tempMatrix(coord.second + (61*i), coord.first + (61*i)) /= currentStationary[coord.first] * proportions[i];
        }
    }

    std::vector<double> dN(truncation, 0.0);
    for(auto coord : nonsynonymous){
        for(int i = 0; i < truncation; i++){
            dN[i] += tempMatrix(coord.first + (i*61), coord.second + (i*61)) * currentStationary[coord.first];
            dN[i] += tempMatrix(coord.second + (i*61), coord.first + (i*61)) * currentStationary[coord.second];
        }
    }

    std::vector<double> dS(truncation, 0.0);
    for(auto coord : synonymous){
        for(int i = 0; i < truncation; i++){
            dS[i] += tempMatrix(coord.first + (i*61), coord.second + (i*61)) * currentStationary[coord.first];
            dS[i] += tempMatrix(coord.second + (i*61), coord.first + (i*61)) * currentStationary[coord.second];
        }
    }

    std::vector<double> dNdSVec(truncation, 0.0);

    for(int i = 0; i < truncation; i++){
        dNdSVec[i] = dN[i]/dS[i];
    }

    return dNdSVec;
}

void SBMatrix::tune(){
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

    if(proportionCount > 0){
        double proportionRate = (double)proportionAcceptCount/(double)proportionCount;

        if ( proportionRate > 0.33 ) {
            proportionAlpha /= (1.0 + ((proportionRate-0.33)/0.67));
        }
        else {
            proportionAlpha *= (2.0 - proportionRate/0.33);
        }

        proportionAcceptCount= 0;
        proportionCount = 0;
    }
}