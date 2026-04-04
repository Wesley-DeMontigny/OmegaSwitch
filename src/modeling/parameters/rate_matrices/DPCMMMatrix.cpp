#include "misc/Math.hpp"
#include "misc/Matrix.hpp"
#include "misc/Probability.hpp"
#include "misc/RandomVariable.hpp"
#include "misc/Settings.hpp"
#include "DPCMMMatrix.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

DPCMMMatrix::DPCMMMatrix(Settings& settings, int nC) : 
                                   currentQMatrix(183, 183, 0.0), oldQMatrix(183, 183, 0.0), currentStationary(61, -1), oldStationary(61, -1), 
                                   kLambda(settings.kLambda), rLambda(settings.rLambda),
                                   stationaryPriorAlpha(61, 2.0), numChar(nC), omegaLambda(settings.omegaLambda), fixedRegimes(settings.fixedRegimes), fixedAssignments(settings.fixCorrectDP && !settings.fixedDPAssignments.empty()), assignments(nC, 0), randomStates(61, 0.0) {


    RandomVariable& rng = RandomVariable::randomVariableInstance();

    currentParams[0] = Probability::Exponential::rv(&rng, kLambda);
    oldParams[0] = currentParams[0];
    currentParamPriors[0] = Probability::Exponential::lnPdf(kLambda, currentParams[0]);
    oldParamPriors[0] = currentParamPriors[0];

    currentParams[1] = Probability::Exponential::rv(&rng, rLambda);
    oldParams[1] = currentParams[1];
    currentParamPriors[1] = Probability::Exponential::lnPdf(rLambda, currentParams[1]);
    oldParamPriors[1] = currentParamPriors[1];

    Probability::Dirichlet::rv(&rng, stationaryPriorAlpha, currentStationary);
    oldStationary = currentStationary;
    currentStationaryPrior = Probability::Dirichlet::lnPdf(stationaryPriorAlpha, currentStationary);
    oldStationaryPrior = currentStationaryPrior;

    if(fixedRegimes > 0){
        currentActiveOmegas = fixedRegimes;
        oldActiveOmegas = fixedRegimes;
    }

    dpAlpha = calculateAlpha(settings.expectedCat, numChar);

    if(fixedAssignments){
        if(settings.fixedDPAssignments.size() != numChar){
            Msg::error("Expected " + std::to_string(numChar) + " fixed DPP assignments, but found " + std::to_string(settings.fixedDPAssignments.size()) + ".");
        }

        int maxAssignment = *std::max_element(settings.fixedDPAssignments.begin(), settings.fixedDPAssignments.end());
        if(maxAssignment < 0){
            Msg::error("Fixed DPP assignments must be non-negative.");
        }

        currentCategories.reserve(maxAssignment + 1);
        for(int i = 0; i <= maxAssignment; i++){
            double omega1 = Probability::Exponential::rv(&rng, omegaLambda);
            double omegaInc1 = Probability::Exponential::rv(&rng, omegaLambda);
            double omegaInc2 = Probability::Exponential::rv(&rng, omegaLambda);
            currentCategories.push_back(Category{{omega1, omegaInc1, omegaInc2}, {}, true});
        }

        for(int i = 0; i < numChar; i++){
            int assignment = settings.fixedDPAssignments[i];
            if(assignment < 0 || assignment >= currentCategories.size()){
                Msg::error("Encountered an invalid fixed DPP assignment.");
            }
            currentCategories[assignment].members.push_back(i);
        }
    }
    else{
        for(int i = 0; i < numChar; i++){
            double randomVal = rng.uniformRv();
            double total = dpAlpha/(i + dpAlpha);

            // If new category
            if(total > randomVal){
                double omega1 = Probability::Exponential::rv(&rng, omegaLambda);
                double omegaInc1 = Probability::Exponential::rv(&rng, omegaLambda);
                double omegaInc2 = Probability::Exponential::rv(&rng, omegaLambda);
                Category newCat = {omega1, omegaInc1, omegaInc2, {i}, true};
                currentCategories.push_back(newCat);
                continue;
            }

            for(Category &c : currentCategories){
                total += c.members.size()/(i+dpAlpha);

                //If old category
                if(total > randomVal){
                    c.members.push_back(i);
                    break;
                }
            }
        }
    }

    oldCategories = currentCategories;

    regenerateAssignments();
    regenerateCatPrior();
    rebuildQMatrix();
    
    oldQMatrix = currentQMatrix.copy();

    std::iota(randomStates.begin(), randomStates.end(), 0);

    dirty();
}

void DPCMMMatrix::rebuildQMatrix(){
    currentQMatrix = Matrix<double>(61*currentActiveOmegas, 61*currentActiveOmegas, 0.0);
    for(const auto& [c1, c2] : MatrixHelper::validPairs){
        for(int i = 0; i < currentActiveOmegas; i++){
            currentQMatrix(c1 + (i*61), c2 + (i*61)) = currentStationary[c2];
            currentQMatrix(c2 + (i*61), c1 + (i*61)) = currentStationary[c1];
        }
    }
    for(const auto& [c1, c2] : MatrixHelper::transitionPairs){
        for(int i = 0; i < currentActiveOmegas; i++){
            currentQMatrix(c1 + (i*61), c2 + (i*61)) *= currentParams[0];
            currentQMatrix(c2 + (i*61), c1 + (i*61)) *= currentParams[0];
        }
    }
}

void DPCMMMatrix::accept() {
    oldParams[0] = currentParams[0];
    oldParamPriors[0] = currentParamPriors[0];
    oldParams[1] = currentParams[1];
    oldParamPriors[1] = currentParamPriors[1];
    oldStationary = currentStationary;
    oldStationaryPrior = currentStationaryPrior;
    oldCatPrior = currentCatPrior;
    oldCategories = currentCategories;
    oldActiveOmegas = currentActiveOmegas;
    oldQMatrix = currentQMatrix.copy();

    if(moveChoice == MatrixMoves::K_MOVE){
        kAcceptCount += 1;
        if(countTuningEvents){
            tuningState->kStats.acceptCount += 1;
        }
    }
    else if(moveChoice == MatrixMoves::STATIONARY_MOVE){
        stationaryAcceptCount += 1;
        if(countTuningEvents){
            tuningState->stationaryStats.acceptCount += 1;
        }
    }
    else if(moveChoice == MatrixMoves::R_MOVE){
        rAcceptCount += 1;
        if(countTuningEvents){
            tuningState->rStats.acceptCount += 1;
        }
    }
    else if(moveChoice == MatrixMoves::OMEGA_MOVE){
        omegaAcceptCount += 1;
        if(countTuningEvents){
            tuningState->omegaStats.acceptCount += 1;
        }
    }

    moveChoice = MatrixMoves::NO_MOVE;
}

void DPCMMMatrix::reject() {
    currentParams[0] = oldParams[0];
    currentParamPriors[0] = oldParamPriors[0];
    currentParams[1] = oldParams[1];
    currentParamPriors[1] = oldParamPriors[1];
    currentStationary = oldStationary;
    currentStationaryPrior = oldStationaryPrior;
    currentCatPrior = oldCatPrior;
    currentCategories = oldCategories;
    currentActiveOmegas = oldActiveOmegas;
    currentQMatrix = oldQMatrix.copy();

    moveChoice = MatrixMoves::NO_MOVE;
}

double DPCMMMatrix::lnPrior() {
    double prior = currentParamPriors[0] + currentStationaryPrior + currentCatPrior;

    if(currentActiveOmegas > 1)
        prior += currentParamPriors[1];

    return prior;
}

void DPCMMMatrix::regenerateCatPrior(){
    int numCats = currentCategories.size();
    
    currentCatPrior = std::log(dpAlpha) * numCats;

    for(Category& c : currentCategories) {
        currentCatPrior += Math::lnFactorial(c.members.size() - 1);
        for(int i = 0; i < currentActiveOmegas; i++){
            currentCatPrior += Probability::Exponential::lnPdf(omegaLambda, c.omegas[i]);
        }
    }
}

double DPCMMMatrix::updateK() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;

    moveChoice = MatrixMoves::K_MOVE;
    kCount += 1;
    if(countTuningEvents){
        tuningState->kStats.count += 1;
    }

    double currentV = currentParams[0];
    double scale = std::exp(tuningState->kDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentParams[0] = newV;
    hastings = std::log(scale);

    this->dirty();

    rebuildQMatrix();
    currentParamPriors[0] = Probability::Exponential::lnPdf(kLambda, currentParams[0]);

    return hastings;
}

double DPCMMMatrix::updateR() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    moveChoice = MatrixMoves::R_MOVE;
    rCount += 1;
    if(countTuningEvents){
        tuningState->rStats.count += 1;
    }

    double currentV = currentParams[1];
    double scale = std::exp(tuningState->rDelta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;

    currentParams[1] = newV;
    hastings = std::log(scale);

    currentParamPriors[1] = Probability::Exponential::lnPdf(rLambda, currentParams[1]);

    return hastings;
}

double DPCMMMatrix::updateStationary() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    this->dirty();
    double hastings = 0.0;

    int numElements = 1;
    moveChoice = MatrixMoves::STATIONARY_MOVE;
    stationaryCount += 1;
    if(countTuningEvents){
        tuningState->stationaryStats.count += 1;
    }

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
        alphaForward[i] = (x[i] * tuningState->stationaryAlpha) + 1.0;
    }
    
    Probability::Dirichlet::rv(&rng, alphaForward, z);

    for(int i = 0; i < z.size(); i++) {
        alphaReverse[i] = (z[i] * tuningState->stationaryAlpha) + 1.0;
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

    rebuildQMatrix();

    return hastings;
}

double DPCMMMatrix::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    this->dirty();
    double hastings = 0.0;

    int randomCategory = (int)(rng.uniformRv() * currentCategories.size());
    currentCategories[randomCategory].dirty = true;

    double randomMove = rng.uniformRv();

    if(randomMove < 0.6 || currentActiveOmegas == 1){ // Scale Move
        moveChoice = MatrixMoves::OMEGA_MOVE;
        omegaCount += 1;
        if(countTuningEvents){
            tuningState->omegaStats.count += 1;
        }

        int randomOmega = (int)(rng.uniformRv() * currentActiveOmegas);
        double currentV = currentCategories[randomCategory].omegas[randomOmega];
        double scale = std::exp(tuningState->omegaDelta * (rng.uniformRv() - 0.5));
        double newV = currentV * scale;

        currentCategories[randomCategory].omegas[randomOmega] = newV;
        hastings = std::log(scale);
    }
    else if(randomMove < 0.9){ // Exchange Move
        moveChoice = MatrixMoves::EXCHANGE_MOVE;

        int randomOmega1 = (int)(currentActiveOmegas * rng.uniformRv());
        int randomOmega2 = 0;

        do {
            randomOmega2 = (int)(currentActiveOmegas * rng.uniformRv());
        }
        while(randomOmega2 == randomOmega1);
        
        double forwardU = Probability::Beta::rv(&rng, 1.0, 5.0);
        double originalO1 = currentCategories[randomCategory].omegas[randomOmega1];

        currentCategories[randomCategory].omegas[randomOmega1] *= (1-forwardU);
        currentCategories[randomCategory].omegas[randomOmega2] += forwardU * originalO1;

        double reverseU = (originalO1 - currentCategories[randomCategory].omegas[randomOmega1]) / currentCategories[randomCategory].omegas[randomOmega2];

        hastings = Probability::Beta::lnPdf(1.0, 5.0, reverseU) - Probability::Beta::lnPdf(1.0, 5.0, forwardU);
        hastings += std::log(originalO1); // Add the Jacobain correction |(1-U, -O1), (U, O1)| = 01
    }
    else { // Re-Index Move
        moveChoice = MatrixMoves::REINDEX_MOVE;

        int randomOmega1 = (int)(currentActiveOmegas * rng.uniformRv());
        int randomOmega2 = 0;

        do {
            randomOmega2 = (int)(currentActiveOmegas * rng.uniformRv());
        }
        while(randomOmega2 == randomOmega1);

        double tempO = currentCategories[randomCategory].omegas[randomOmega1];

        currentCategories[randomCategory].omegas[randomOmega1] = currentCategories[randomCategory].omegas[randomOmega2];
        currentCategories[randomCategory].omegas[randomOmega2] = tempO;
    }

    regenerateCatPrior();

    return hastings;
}


double DPCMMMatrix::updateActiveOmegas() {
    if(isRegimeCountFixed()){
        return -1 * INFINITY;
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hastings = 0.0;
    this->dirty();

    double randomMove = rng.uniformRv();

    if(randomMove > 0.5){
        moveChoice = MatrixMoves::SPLIT_MERGE_MOVE;
        double splitAlpha = 5.0;
        double total = MatrixHelper::possibleSplit3[currentActiveOmegas-1] + MatrixHelper::possibleMerge3[currentActiveOmegas - 1];
        double splitProbs = (double)(MatrixHelper::possibleSplit3[currentActiveOmegas-1]) / total;
        double randomMove2 = rng.uniformRv();

        if(randomMove2 < splitProbs){ // Perform a split
            double forwardProb = -std::log(total);
            double nextTotal = MatrixHelper::possibleSplit3[currentActiveOmegas] + MatrixHelper::possibleMerge3[currentActiveOmegas];
            double reverseProb = -std::log(nextTotal);
            hastings += reverseProb - forwardProb;

            int randomSplit = (int)(rng.uniformRv() * currentActiveOmegas);

            for (int i = 0; i < currentCategories.size(); i++) {
                double u = Probability::Beta::rv(&rng, splitAlpha, splitAlpha);
                
                // If we select a middle index to split we need to shift everything over
                for(int j = currentActiveOmegas; j > randomSplit + 1; j--){
                    currentCategories[i].omegas[j] = currentCategories[i].omegas[j-1];
                }

                double tempO = currentCategories[i].omegas[randomSplit];

                currentCategories[i].omegas[randomSplit] = tempO * u;
                currentCategories[i].omegas[randomSplit + 1] = tempO * (1.0 - u);
                currentCategories[i].dirty = true;

                hastings += std::log(tempO);
                hastings -= Probability::Beta::lnPdf(splitAlpha, splitAlpha, u);
            }

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
            double nextTotal = MatrixHelper::possibleSplit3[currentActiveOmegas - 2] + MatrixHelper::possibleMerge3[currentActiveOmegas - 2];
            double reverseProb = -std::log(nextTotal);
            hastings += reverseProb - forwardProb;

            int randomMerge = (int)(rng.uniformRv() * (currentActiveOmegas - 1));
            
            for(int i = 0; i < currentCategories.size(); i++){
                double o1 = currentCategories[i].omegas[randomMerge];
                double o2 = currentCategories[i].omegas[randomMerge + 1];
                double total = o1 + o2;
                double u = o1 / total;

                double sum = currentCategories[i].omegas[randomMerge] + currentCategories[i].omegas[randomMerge + 1];
                currentCategories[i].omegas[randomMerge] = sum;
                currentCategories[i].dirty = true;

                // We need to shift elements down depending on the random merge
                for(int j = randomMerge + 1; j < currentActiveOmegas-1; j++){
                    currentCategories[i].omegas[j] = currentCategories[i].omegas[j + 1];
                }

                hastings += Probability::Beta::lnPdf(splitAlpha, splitAlpha, u);
                hastings -= std::log(sum);
            }

            currentActiveOmegas--;

            if(currentActiveOmegas == 1){
                hastings += Probability::Exponential::lnPdf(rLambda, currentParams[1]);
            }
        }
    }
    // We don't need a hastings correction for the chosen index. If we have n omegas, there will be n choices for insertion after I delete and have n-1 omegas left
    else{ // Pure Birth-Death Moves
        moveChoice = MatrixMoves::BIRTH_DEATH_MOVE;
        
        if((rng.uniformRv() > 0.5 || currentActiveOmegas == 1) && currentActiveOmegas != 3){ // Pure Birth
            int randomBirth = (int)(rng.uniformRv() * (currentActiveOmegas + 1));

            for (int i = 0; i < currentCategories.size(); i++) {
                for(int j = currentActiveOmegas; j > randomBirth; j--){
                    currentCategories[i].omegas[j] = currentCategories[i].omegas[j-1];
                }

                double birth = Probability::Exponential::rv(&rng, omegaLambda);

                currentCategories[i].omegas[randomBirth] = birth;
                currentCategories[i].dirty = true;

                hastings -= Probability::Exponential::lnPdf(omegaLambda, birth);
            }

            currentActiveOmegas++;

            if(currentActiveOmegas == 2){
                double newR = Probability::Exponential::rv(&rng, rLambda);
                currentParamPriors[1] = Probability::Exponential::lnPdf(rLambda, newR);
                hastings -= currentParamPriors[1];

                currentParams[1] = newR;

                hastings -= std::log(2.0); // Correct for forced birth
            }
            else if(currentActiveOmegas == 3){
                hastings += std::log(2.0); // Correct for forced death in the future
            }
        }
        else{ // Pure Death
            int randomDeath = (int)(rng.uniformRv() * currentActiveOmegas);

            for(int i = 0; i < currentCategories.size(); i++){
                hastings += Probability::Exponential::lnPdf(omegaLambda, currentCategories[i].omegas[randomDeath]);

                for(int j = randomDeath; j < currentActiveOmegas-1; j++){
                    currentCategories[i].omegas[j] = currentCategories[i].omegas[j + 1];
                }

                currentCategories[i].dirty = true;
            }

            currentActiveOmegas--;

            if(currentActiveOmegas == 1){
                hastings += Probability::Exponential::lnPdf(rLambda, currentParams[1]);
                hastings += std::log(2.0); // Correct for forced birth in the future
            }
            else if(currentActiveOmegas == 2){
                hastings -= std::log(2.0); // Correct for forced death
            }
        }
    }
    

    rebuildQMatrix();
    regenerateCatPrior();

    return hastings;
}


double DPCMMMatrix::expectedCategories(double a, int members) const {
    return a * std::log(1 + (members/a));
}

// From John's code
double DPCMMMatrix::calculateAlpha(double expectedCat, int members) const {

    if (expectedCat > members)
        Msg::error("The expected number of tables cannot be larger than the number of patrons (" + std::to_string(members) + "<" + std::to_string(expectedCat) + ")");
     if (expectedCat < 1.0)
        Msg::error("The expected number of tables cannot be less than one");
       
    double a = 0.000001;
    double ea = expectedCategories(a, members);
    bool goUp;
    if (ea < expectedCat)
        goUp = true;
    else
        goUp = false;
    double increment = 0.1;
    
    //While we are not within some error tolerance, move increment the estimate to get closer and closer
    while (fabs(ea - expectedCat) > 0.000001) {
        if (ea < expectedCat && goUp == true){
            a += increment;
        }
        else if (ea > expectedCat && goUp == false){
            a -= increment;
        }
        else if (ea < expectedCat && goUp == false){
            increment /= 2.0;
            goUp = true;
            a += increment;
        }
        else{
            increment /= 2.0;
            goUp = false;
            a -= increment;
        }
        ea = expectedCategories(a, members);
    }
    return a;
}

void DPCMMMatrix::assignMember(int member, int category){
    currentCategories[category].members.push_back(member);
}

void DPCMMMatrix::regenerateAssignments(){
    int numCats = currentCategories.size();
    int checkSum = 0;
    for(int i = 0; i < numCats; i++){
        for(int m : currentCategories[i].members){
            assignments[m] = i;
        }
        checkSum += currentCategories[i].members.size();
    }
    if(checkSum != numChar)
        Msg::error("Assignment book keeping failed!");
}

void DPCMMMatrix::removeCategory(int index){
    if(currentCategories[index].members.size() != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    currentCategories.erase(currentCategories.begin() + index);
}

void DPCMMMatrix::addCategory(const std::array<double, 3>& omegas){
    Category newCat = {omegas, {}, true};
    currentCategories.push_back(newCat);
}

std::optional<int> DPCMMMatrix::unassignMember(int member){
    for(int i = 0; i < currentCategories.size(); i++){
        Category& c = currentCategories[i];
        for(int j = 0; j < c.members.size(); j++){
            if(c.members[j] == member){
                c.members.erase(c.members.begin() + j);
                if(c.members.size() == 0){
                    removeCategory(i);
                    return i;
                }
                return std::nullopt;
            }
        }
    }
    return std::nullopt;
}

void DPCMMMatrix::popCategories(int n){
    int size = currentCategories.size();
    for(int i = 0; i < n; i++){
        currentCategories.pop_back();
    }
}

Matrix<double> DPCMMMatrix::Q(int c) const {
    return Q(currentCategories[c].omegas);
}

Matrix<double> DPCMMMatrix::Q(const std::array<double, 3>& omegas) const {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    int stateSpace = currentQMatrix.dim1();

    for(const auto& [c1, c2] : MatrixHelper::nonsynonymousPairs){
        double total = 0.0;
        for(int i = 0; i < currentActiveOmegas; i++){
            total += omegas[i];
            returnMatrix(c1 + (i*61), c2 + (i*61)) *= total;
            returnMatrix(c2 + (i*61), c1 + (i*61)) *= total;
        }
    }

    double stateSpaceScaler = 1.0 / (double)(currentActiveOmegas);
    double scaler= 0.0;
    for(int i = 0; i < stateSpace; i++){
        double total = 0.0;
        for(int j = 0; j < stateSpace; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1.0;
        scaler += returnMatrix(i, i) * currentStationary[i % 61] * stateSpaceScaler;
    }

    scaler = -1.0 / scaler;
    for (int i = 0; i < stateSpace; i++)
        for (int j = 0; j < stateSpace; j++)
            returnMatrix(i, j) *= scaler;
    
    if(stateSpace == 183){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = currentParams[1];
            returnMatrix(i + 61, i) = currentParams[1];
            returnMatrix(i, i + 122) = currentParams[1];
            returnMatrix(i + 122, i) = currentParams[1];
            returnMatrix(i + 61, i + 122) = currentParams[1];
            returnMatrix(i + 122, i + 61) = currentParams[1];

            returnMatrix(i, i) -= 2.0*currentParams[1];
            returnMatrix(i + 61, i + 61) -= 2.0*currentParams[1];
            returnMatrix(i + 122, i + 122) -= 2.0*currentParams[1];
        }
    }
    else if(stateSpace == 122){
        for(int i = 0; i < 61; i++){
            returnMatrix(i, i + 61) = currentParams[1];
            returnMatrix(i + 61, i) = currentParams[1];

            returnMatrix(i, i) -= currentParams[1];
            returnMatrix(i + 61, i + 61) -= currentParams[1];
        } 
    }

    return returnMatrix;
}

std::vector<double> DPCMMMatrix::getStationary(){
    std::vector<double> returnStationary;

    double factor = 1.0/(double)(currentActiveOmegas);

    for(int i = 0; i < currentActiveOmegas; i++){
        for(double v : currentStationary){
            returnStationary.push_back(v * factor);
        }
    }
    
    return returnStationary;
}

std::array<double, 3> DPCMMMatrix::dNdS(int c) const {
    Matrix<double> tempMatrix(183, 183, 0.0);
    std::array<double, 3> dNdS = {0,0,0}; // We start out with dN and then divide by dS
    std::array<double, 3> dS = {0,0,0};

    for(const auto& [c1, c2] : MatrixHelper::validPairs){
        for(int i = 0; i < 3; i++){
            tempMatrix(c1 + (i*61), c2 + (i*61)) = 1.0;
            tempMatrix(c2 + (i*61), c1 + (i*61)) = 1.0;
        }
    }
    for(const auto& [c1, c2] : MatrixHelper::nonsynonymousPairs){
        double total = 0.0;
        for(int i = 0; i < 3; i++){
            total += currentCategories[c].omegas[i];
            tempMatrix(c1 + (i*61), c2 + (i*61)) *= total;
            tempMatrix(c2 + (i*61), c1 + (i*61)) *= total;
        }
    }
    for(const auto& [c1, c2] : MatrixHelper::transitionPairs){
        for(int i = 0; i < 3; i++){
            tempMatrix(c1 + (i*61), c2 + (i*61)) *= currentParams[0];
            tempMatrix(c2 + (i*61), c1 + (i*61)) *= currentParams[0];
        }
    }

    for(const auto& [c1, c2] : MatrixHelper::nonsynonymousPairs){
        for(int i = 0; i < 3; i++){
            dNdS[i] += tempMatrix(c1 + (i*61), c2 + (i*61)) * currentStationary[c1];
            dNdS[i] += tempMatrix(c2 + (i*61), c1 + (i*61)) * currentStationary[c2];
        }
    }

    for(const auto& [c1, c2] : MatrixHelper::synonymousPairs){
        for(int i = 0; i < 3; i++){
            dS[i] += tempMatrix(c1 + (i*61), c2 + (i*61)) * currentStationary[c1];
            dS[i] += tempMatrix(c2 + (i*61), c1 + (i*61)) * currentStationary[c2];
        }
    }

    for(int i = 0; i < 3; i++){
        dNdS[i] /= dS[i];
    }

    return dNdS;
}

double DPCMMMatrix::getOmegaRate() const {
    return (double)omegaAcceptCount/(double)omegaCount;
}

double DPCMMMatrix::getStationaryRate() const {
    return (double)stationaryAcceptCount/(double)stationaryCount;
}

double DPCMMMatrix::getRRate() const {
    return (double)rAcceptCount/(double)rCount;
}

double DPCMMMatrix::getKRate() const {
    return (double)kAcceptCount/(double)kCount;
}

void DPCMMMatrix::tune(){
    if(tuningState->kStats.count > 0){
        double kRate = (double)tuningState->kStats.acceptCount / (double)tuningState->kStats.count;

        if(kRate > 0.33){
            tuningState->kDelta *= (1.0 + ((kRate-0.33)/0.67));
        }
        else{
            tuningState->kDelta /= (2.0 - kRate/0.33);
        }
        tuningState->kStats.acceptCount = 0;
        tuningState->kStats.count = 0;
    }

    if(tuningState->rStats.count > 0){
        double rRate = (double)tuningState->rStats.acceptCount / (double)tuningState->rStats.count;

        if(rRate > 0.33){
            tuningState->rDelta *= (1.0 + ((rRate-0.33)/0.67));
        }
        else{
            tuningState->rDelta /= (2.0 - rRate/0.33);
        }
        tuningState->rStats.acceptCount = 0;
        tuningState->rStats.count = 0;
    }

    if(tuningState->stationaryStats.count > 0){
        double stationaryRate = (double)tuningState->stationaryStats.acceptCount / (double)tuningState->stationaryStats.count;

        if(stationaryRate > 0.33){
            tuningState->stationaryAlpha /= (1.0 + ((stationaryRate-0.33)/0.67));
        }
        else{
            tuningState->stationaryAlpha *= (2.0 - stationaryRate/0.33);
        }

        tuningState->stationaryStats.acceptCount = 0;
        tuningState->stationaryStats.count = 0;
    }

    if(tuningState->omegaStats.count > 0){
        double omegaRate = (double)tuningState->omegaStats.acceptCount / (double)tuningState->omegaStats.count;

        if(omegaRate > 0.33){
            tuningState->omegaDelta *= (1.0 + ((omegaRate-0.33)/0.67));
        }
        else{
            tuningState->omegaDelta /= (2.0 - omegaRate/0.33);
        }
        tuningState->omegaStats.acceptCount = 0;
        tuningState->omegaStats.count = 0;
    }
}
