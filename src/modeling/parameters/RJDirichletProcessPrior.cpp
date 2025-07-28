#include "RJDirichletProcessPrior.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/RJDPPModel.hpp" // Annoying circular dependency... It is what is is for now...
#include "modeling/parameters/RJDPPMatrix.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/model/ConditionalLikelihood.hpp"
#include "core/Probability.hpp"
#include "core/Msg.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <chrono>
#include "core/Math.hpp"
#include <iostream>
#include <algorithm>

RJDirichletProcessPrior::RJDirichletProcessPrior(RJDPPMatrix* matrix, int size, Settings s) : 
                                             alpha(0), omegaLambda(s.omegaLambda),
                                             numMembers(size), currentLnPrior(0.0),
                                             model(nullptr), omegaDelta(0.5), assignments(size, -1),
                                             omegaAcceptCount(0), omegaCount(0), moveChoice(-1),
                                             numGibbs(s.numGibbs), rateMatrix(matrix) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    alpha = calculateAlpha(s.expectedCat, numMembers);

    for(int i = 0; i < size; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            double omega1 = Probability::Exponential::rv(&rng, omegaLambda);
            double omegaInc1 = Probability::Exponential::rv(&rng, omegaLambda);
            double omegaInc2 = Probability::Exponential::rv(&rng, omegaLambda);
            double omega2 = omega1 + omegaInc1;
            double omega3 = omega2 + omegaInc2;
            RJCategory newCat = {omega1, omega2, omega3,
                               omegaInc1, omegaInc2,
                               1, {i}, true};
            currentCategories.push_back(newCat);
            continue;
        }

        for(RJCategory &c : currentCategories){
            total += c.size/(i+alpha);

            //If old category
            if(total > randomVal){
                c.size++;
                c.members.push_back(i);
                break;
            }
        }
    }

    int numCats = currentCategories.size();
    for(int i = 0; i < numCats; i++)
        for(int m : currentCategories[i].members)
            assignments[m] = i;
    
    this->dirty();

    regeneratePrior();

    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;
}

RJDirichletProcessPrior::~RJDirichletProcessPrior() {}

void RJDirichletProcessPrior::registerModel(RJDPPModel* m) { model = m; }

double RJDirichletProcessPrior::expectedCategories(double a, int members){
    return a * std::log(1 + (members/a));
}

// From John's code
double RJDirichletProcessPrior::calculateAlpha(double expectedCat, int members) {

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

void RJDirichletProcessPrior::regeneratePrior(){
    int numCats = currentCategories.size();
    
    currentLnPrior = std::log(alpha) * numCats;

    for(RJCategory& c : currentCategories) {
        currentLnPrior += Math::lnFactorial(c.size - 1);
        if(currentActiveOmegas == 3){
            currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega1);
            currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omegaIncrement1);
            currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omegaIncrement2);  
        }
        else if(currentActiveOmegas == 2){
            currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega1);
            currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omegaIncrement1);
        }
        else if(currentActiveOmegas == 1){
            currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega1);
        }
    }
}

void RJDirichletProcessPrior::removeCategory(int index){
    if(currentCategories[index].size != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    currentCategories.erase(currentCategories.begin() + index);
}

void RJDirichletProcessPrior::addCategory(double omega1, double omegaInc1, double omegaInc2){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    RJCategory newCat = {omega1, omega1 + omegaInc1, omega1 + omegaInc1 + omegaInc2, omegaInc1, omegaInc2, 0, {}, true};
    currentCategories.push_back(newCat);
}

int RJDirichletProcessPrior::unassignMember(int member){
    for(int i = 0; i < currentCategories.size(); i++){
        RJCategory& c = currentCategories[i];
        for(int j = 0; j < c.size; j++){
            if(c.members[j] == member){
                c.size--;
                c.members.erase(c.members.begin() + j);
                if(c.size == 0){
                    removeCategory(i);
                    return i;
                }
                return -1;
            }
        }
    }
    return -2;
}

void RJDirichletProcessPrior::assignMember(int member, int category){
    currentCategories[category].size++;
    currentCategories[category].members.push_back(member);
}

void RJDirichletProcessPrior::accept() {
    if(moveChoice == 0){
        omegaAcceptCount += 1;
    }

    moveChoice = -1;

    for(int i = 0; i < currentCategories.size(); i++){
        currentCategories[i].members.shrink_to_fit();
        currentCategories[i].dirty = false;
    }

    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;
    oldActiveOmegas = currentActiveOmegas;
}

void RJDirichletProcessPrior::reject() {
    currentCategories = oldCategories;
    currentLnPrior = oldLnPrior;
    currentActiveOmegas = oldActiveOmegas;

    moveChoice = -1;
}

double RJDirichletProcessPrior::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    this->dirty();
    double hastings = 0.0;

    int randomCategory = (int)(rng.uniformRv() * currentCategories.size());
    currentCategories[randomCategory].dirty = true;

    moveChoice = 0;
    omegaCount += 1;

    if(currentActiveOmegas == 3){
        int randomOmega = (int)(rng.uniformRv() * 3);
        if(randomOmega == 0){
            double currentV1 = currentCategories[randomCategory].omega1;
            double scale1 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV1 = currentV1 * scale1;

            currentCategories[randomCategory].omega1 = newV1;
            currentCategories[randomCategory].omega2 = currentCategories[randomCategory].omegaIncrement1 + newV1;
            currentCategories[randomCategory].omega3 = currentCategories[randomCategory].omegaIncrement1 + currentCategories[randomCategory].omegaIncrement2 + newV1;
            hastings = std::log(scale1);
        }
        else if(randomOmega == 1) {
            double currentV2 = currentCategories[randomCategory].omegaIncrement1;
            double scale2 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV2 = currentV2 * scale2;

            currentCategories[randomCategory].omegaIncrement1 = newV2;
            currentCategories[randomCategory].omega2 = currentCategories[randomCategory].omega1 + newV2;
            currentCategories[randomCategory].omega3 = currentCategories[randomCategory].omega1 + currentCategories[randomCategory].omegaIncrement2 + newV2;
            hastings = std::log(scale2);
        }
        else if(randomOmega == 2) {
            double currentV3 = currentCategories[randomCategory].omegaIncrement2;
            double scale3 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV3 = currentV3 * scale3;

            currentCategories[randomCategory].omegaIncrement2 = newV3;
            currentCategories[randomCategory].omega3 = currentCategories[randomCategory].omega1 + currentCategories[randomCategory].omegaIncrement1 + newV3;
            hastings = std::log(scale3);
        }
    }
    else if(currentActiveOmegas == 2){
        int randomOmega = (int)(rng.uniformRv() * 2);
        if(randomOmega == 0){
            double currentV1 = currentCategories[randomCategory].omega1;
            double scale1 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV1 = currentV1 * scale1;

            currentCategories[randomCategory].omega1 = newV1;
            currentCategories[randomCategory].omega2 = currentCategories[randomCategory].omegaIncrement1 + newV1;
            hastings = std::log(scale1);
        }
        else if(randomOmega == 1) {
            double currentV2 = currentCategories[randomCategory].omegaIncrement1;
            double scale2 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
            double newV2 = currentV2 * scale2;

            currentCategories[randomCategory].omegaIncrement1 = newV2;
            currentCategories[randomCategory].omega2 = currentCategories[randomCategory].omega1 + newV2;
            hastings = std::log(scale2);
        }
    }
    else if(currentActiveOmegas == 1){
        double currentV1 = currentCategories[randomCategory].omega1;
        double scale1 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV1 = currentV1 * scale1;

        currentCategories[randomCategory].omega1 = newV1;
        hastings = std::log(scale1);
    }

    regeneratePrior();

    return hastings;
}

double RJDirichletProcessPrior::updateActiveOmegas() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double alpha = 15.0;
    double hastings = 0.0;
    this->dirty();

    if (currentActiveOmegas == 1) { // 1 2 split
        hastings = std::log(2.0);
        currentActiveOmegas = 2;

        double u = Probability::Beta::rv(&rng, alpha, alpha);
        for (int i = 0; i < currentCategories.size(); i++) {
            double tempO = currentCategories[i].omega1;
            currentCategories[i].omega1 = tempO * u;
            currentCategories[i].omegaIncrement1 = tempO * (1.0 - u);
            currentCategories[i].omega2 = currentCategories[i].omegaIncrement1 + currentCategories[i].omega1;
            currentCategories[i].dirty = true;
            hastings += std::log(tempO);
        }

        hastings += Probability::Beta::lnPdf(alpha, alpha, u);
        
        double newR = Probability::Exponential::rv(&rng, rateMatrix->getRLambda());
        hastings += Probability::Exponential::lnPdf(rateMatrix->getRLambda(), newR);

        rateMatrix->setR(newR);
    }
    else if (currentActiveOmegas == 2) { // 2 1 merge or 2 3 split
        hastings = std::log(0.5);

        if (rng.uniformRv() > 0.5) { // 2 3 split
            currentActiveOmegas = 3;

            double u = Probability::Beta::rv(&rng, alpha, alpha);
            for (int i = 0; i < currentCategories.size(); i++) {
                double tempO = currentCategories[i].omegaIncrement1;
                currentCategories[i].omegaIncrement1 = tempO * u;
                currentCategories[i].omegaIncrement2 = tempO * (1.0 - u);
                currentCategories[i].omega2 = currentCategories[i].omegaIncrement1 + currentCategories[i].omega1;
                currentCategories[i].omega2 = currentCategories[i].omegaIncrement2 + currentCategories[i].omega2;
                currentCategories[i].dirty = true;
                hastings += std::log(tempO);
            }

            hastings += Probability::Beta::lnPdf(alpha, alpha, u);
        } else { // 2 1 merge
            currentActiveOmegas = 1;
            double o1 = currentCategories[0].omega1;
            double o2 = currentCategories[0].omegaIncrement1;
            double total = o1 + o2;
            double u = o1 / total;
            hastings -= Probability::Beta::lnPdf(alpha, alpha, u); // It is identical for all categories

            for(int i = 0; i < currentCategories.size(); i++){
                double sum = currentCategories[i].omega1 + currentCategories[i].omegaIncrement1;
                currentCategories[i].omega1 = sum;
                currentCategories[i].dirty = true;
                hastings -= std::log(sum);
            }
        }
    }
    else if (currentActiveOmegas == 3) { // 3 2 merge
        currentActiveOmegas = 2;
        hastings = std::log(0.5);

        double o2 = currentCategories[0].omegaIncrement1;
        double o3 = currentCategories[0].omegaIncrement2;
        double total = o2 + o3;
        double u = o2 / total;
        hastings -= Probability::Beta::lnPdf(alpha, alpha, u);

        for(int i = 0; i < currentCategories.size(); i++){
            double sum = currentCategories[i].omegaIncrement1 + currentCategories[i].omegaIncrement2;
            currentCategories[i].omegaIncrement1 = sum;
            currentCategories[i].omega2 = currentCategories[i].omegaIncrement1 + currentCategories[i].omega1;
            currentCategories[i].dirty = true;
            hastings -= std::log(sum);
        }
    }

    return hastings;
}

double RJDirichletProcessPrior::updateDPP(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    this->dirty();

    int numAux = 5;

    for(int iter = 0; iter < numGibbs; iter++) {
        #if TIME_PROFILE==1
        std::chrono::steady_clock::time_point initTime = std::chrono::steady_clock::now();
        #endif
        int n = (int)(rng.uniformRv() * numMembers);

        int deleted = unassignMember(n); // This will also delete the group if empty
        if(deleted >= 0){
            model->getTransitionProbability()->deleteQ(deleted);
        }

        std::vector<double> conditionalL;
        int numCats = currentCategories.size();

        tf::Taskflow taskflow;

        for(int i = 0; i < numCats; i++){
            conditionalL.push_back(0.0);

            taskflow.emplace([this, &conditionalL, i, n](){
                double likelihood = model->testCategory(n, i, false);
                conditionalL[i] = likelihood + std::log(currentCategories[i].size);
            });
        }

        std::vector<double> omega1Vec;
        std::vector<double> omegaInc1Vec;
        std::vector<double> omegaInc2Vec;

        double alphaSplit = std::log(alpha/numAux);

        #if TIME_PROFILE==1
        std::chrono::steady_clock::time_point preAllocation = std::chrono::steady_clock::now();
        #endif
        model->getTransitionProbability()->allocateQ(numCats + numAux);
        #if TIME_PROFILE==1
        std::chrono::steady_clock::time_point postAllocation = std::chrono::steady_clock::now();
        std::cout << "Allocation of Q matrices completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(postAllocation - preAllocation).count() << "[milliseconds]" << std::endl;
        #endif
        for(int i = 0; i < numAux; i++){
            conditionalL.push_back(0.0);
            double newOmega1 = Probability::Exponential::rv(&rng, omegaLambda);
            double newOmegaInc1 = Probability::Exponential::rv(&rng, omegaLambda);
            double newOmegaInc2 = Probability::Exponential::rv(&rng, omegaLambda);

            addCategory(newOmega1, newOmegaInc1, newOmegaInc2);
            omega1Vec.push_back(newOmega1);
            omegaInc1Vec.push_back(newOmegaInc1);
            omegaInc2Vec.push_back(newOmegaInc2);

            taskflow.emplace([this, &conditionalL, i, n, numCats, alphaSplit](){
                double likelihood = model->testCategory(n, numCats+i, true);
                conditionalL[numCats + i] = likelihood + alphaSplit;
            });
        }

        executor.run(taskflow).wait();

        for(int i = 0; i < numAux; i++)
            currentCategories.pop_back();

        //Do some adjustments here to get relative probabilities
        double maxL = *std::max_element(conditionalL.begin(), conditionalL.end());
        double total = 0.0;
        for(double& d : conditionalL){
            d -= maxL;
            d = std::exp(d);
            total += d;
        }

        double categoryDraw = total * rng.uniformRv();

        total = 0.0;
        bool assigned = false;
        for(int i = 0; i < conditionalL.size(); i++){
            total += conditionalL[i];
            if(total > categoryDraw){
                if(i < numCats) { //It already exists
                    assignMember(n, i);
                }
                else {
                    addCategory(omega1Vec[i - numCats], omegaInc1Vec[i - numCats], omegaInc2Vec[i - numCats]);
                    assignMember(n, numCats);
                    model->regenerateTransitionProbs(n, numCats);
                }
                break;
            }
        }

        int bufferDiff = model->getTransitionProbability()->getNumMatrices() - currentCategories.size();
        if(bufferDiff > 30){
            model->getTransitionProbability()->deleteNQ(bufferDiff - 30);
        }

        #if TIME_PROFILE==1
        std::chrono::steady_clock::time_point finalTime = std::chrono::steady_clock::now();
        std::cout << "DPP Iteration was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(finalTime - initTime).count() << "[milliseconds]" << std::endl;
        #endif
    }

    regeneratePrior();

    currentCategories.shrink_to_fit();

    for(int i = 0; i < numMembers; i++)
        assignments[i] = -1;

    for(int i = 0; i < currentCategories.size(); i++){
        for(int m : currentCategories[i].members){
            if(assignments[m] != -1)
                Msg::error("Duplicate Assignments!");
            else
            assignments[m] = i;
        }
    }

    for(int i = 0; i < numMembers; i++){
        if(assignments[i] == -1){
            std::cout << i << " was not assigned to any category." << std::endl;
            Msg::error("Failed to assign!");
        }
    }

    return INFINITY;
}

void RJDirichletProcessPrior::tune() {
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