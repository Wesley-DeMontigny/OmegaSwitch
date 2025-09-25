#include "DirichletProcessPrior.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/DPPModel.hpp" // Annoying circular dependency... It is what is is for now...
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

DirichletProcessPrior::DirichletProcessPrior(int size, Settings s) : 
                                             alpha(0), omegaLambda(s.omegaLambda),
                                             numMembers(size), currentLnPrior(0.0),
                                             model(nullptr), omegaDelta(0.5), assignments(size, -1),
                                             omegaAcceptCount(0), omegaCount(0), moveChoice(-1),
                                             numGibbs(s.numGibbs) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    alpha = calculateAlpha(s.expectedCat, numMembers);

    for(int i = 0; i < size; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            double omega1 = Probability::Exponential::rv(&rng, omegaLambda);
            double omegaI = Probability::Exponential::rv(&rng, omegaLambda);
            double omega2 = omega1 + omegaI;
            Category newCat = {omega1,
                               omega2,
                               omegaI,
                               1, {i}, true};
            currentCategories.push_back(newCat);
            continue;
        }

        for(Category &c : currentCategories){
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

DirichletProcessPrior::~DirichletProcessPrior() {}

void DirichletProcessPrior::registerModel(DPPModel* m) { model = m; }

double DirichletProcessPrior::expectedCategories(double a, int members){
    return a * std::log(1 + (members/a));
}

// From John's code
double DirichletProcessPrior::calculateAlpha(double expectedCat, int members) {

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

void DirichletProcessPrior::regeneratePrior(){
    int numCats = currentCategories.size();
    
    currentLnPrior = std::log(alpha) * numCats;

    for(Category& c : currentCategories) {
        currentLnPrior += Math::lnFactorial(c.size - 1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omegaIncrement);
    }
}

void DirichletProcessPrior::removeCategory(int index){
    if(currentCategories[index].size != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    currentCategories.erase(currentCategories.begin() + index);
}

void DirichletProcessPrior::addCategory(double omega1, double omegaI){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Category newCat = {omega1, omega1 + omegaI, omegaI, 0, {}, true};
    currentCategories.push_back(newCat);
}

int DirichletProcessPrior::unassignMember(int member){
    for(int i = 0; i < currentCategories.size(); i++){
        Category& c = currentCategories[i];
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

void DirichletProcessPrior::assignMember(int member, int category){
    currentCategories[category].size++;
    currentCategories[category].members.push_back(member);
}

void DirichletProcessPrior::accept() {
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
}

void DirichletProcessPrior::reject() {
    currentCategories = oldCategories;
    currentLnPrior = oldLnPrior;

    moveChoice = -1;
}

double DirichletProcessPrior::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    this->dirty();
    double hastings = 0.0;

    int randomCategory = (int)(rng.uniformRv() * currentCategories.size());
    currentCategories[randomCategory].dirty = true;

    int randomOmega = (int)(rng.uniformRv() * 2);

    moveChoice = 0;
    omegaCount += 1;

    if(randomOmega == 0){
        double currentV1 = currentCategories[randomCategory].omega1;
        double scale1 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV1 = currentV1 * scale1;

        currentCategories[randomCategory].omega1 = newV1;
        hastings = std::log(scale1);
    }
    else {
        double currentV2 = currentCategories[randomCategory].omegaIncrement;
        double scale2 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV2 = currentV2 * scale2;

        currentCategories[randomCategory].omegaIncrement = newV2;
        currentCategories[randomCategory].omega2 = currentCategories[randomCategory].omega1 + newV2;
        hastings = std::log(scale2);
    }

    regeneratePrior();

    return hastings;
}

double DirichletProcessPrior::updateDPP(){
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
        std::vector<double> omegaIVec;

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
            double newOmegaI = Probability::Exponential::rv(&rng, omegaLambda);

            addCategory(newOmega1, newOmegaI);
            omega1Vec.push_back(newOmega1);
            omegaIVec.push_back(newOmegaI);

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
                    assigned = true;
                }
                else {
                    addCategory(omega1Vec[i - numCats], omegaIVec[i - numCats]);
                    assignMember(n, numCats);
                    model->regenerateTransitionProbs(n, numCats);
                    assigned = true;
                }
                break;
            }
        }

        if(assigned == false){
            std::cout << n << " was not assigned to any category.\n";
            std::cout << "Log Likelihoods: \n";
            for(double& a : conditionalL){
                std::cout << a << "\n";
            }
            std::cout << std::flush;
            Msg::error("Failed to assign!");
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

    return INFINITY;
}

void DirichletProcessPrior::tune() {
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