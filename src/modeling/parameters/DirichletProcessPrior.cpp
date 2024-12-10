#include "DirichletProcessPrior.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp" // Annoying circular dependency... It is what is is for now...
#include "modeling/model/TransitionProbability.hpp"
#include "core/Probability.hpp"
#include "core/Msg.hpp"
#include <cmath>
#include "core/Math.hpp"
#include <iostream>
#include <algorithm>

DirichletProcessPrior::DirichletProcessPrior(int size, double a, double oL, double rL, int numGibbs) : 
                                             alpha(a), omegaLambda(oL), numMembers(size), currentLnPrior(0.0), 
                                             oldLnPrior(0.0), numGibbsUpdate(numGibbs), model(nullptr),
                                             moveChoice(-1), omegaCount(0), omegaAcceptCount(0), omegaDelta(std::log(2)),
                                             betaCount(0), betaAcceptCount(0), betaDelta(1.0), rLambda(rL) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //The expected category number can be roughly computed as n = a * ln(1 + c/a)
    for(int i = 0; i < size; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            Category newCat = {Probability::Exponential::rv(&rng, omegaLambda),
                               Probability::Exponential::rv(&rng, omegaLambda),
                               Probability::Exponential::rv(&rng, rLambda),
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

    for(int i = 0; i < numMembers; i++)
        assignments.push_back(-1);

    denominator = Math::lnGamma(numMembers - alpha);

    int numCats = currentCategories.size();
    for(int i = 0; i < numCats; i++)
        for(int m : currentCategories[i].members)
            assignments[m] = i;
    
    regeneratePrior();

    oldAssignments = assignments;
    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;
}

DirichletProcessPrior::~DirichletProcessPrior() {
    
}

void DirichletProcessPrior::regeneratePrior(){
    int numCats = currentCategories.size();
    
    currentLnPrior = std::log(alpha) * numCats;

    for (int i = 0; i < numCats; ++i) {
        currentLnPrior += Math::lnFactorial(currentCategories[i].size - 1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, currentCategories[i].omega);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, currentCategories[i].beta);
        currentLnPrior += Probability::Exponential::lnPdf(rLambda, currentCategories[i].r);
    }

    currentLnPrior -= denominator;
}

void DirichletProcessPrior::removeCategory(int index){
    if(currentCategories[index].size != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    currentCategories.erase(currentCategories.begin() + index);
}

void DirichletProcessPrior::addCategory(double omega1, double omega2, double r){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Category newCat = {omega1, omega2, r, 0, {}, true};
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

int DirichletProcessPrior::unassignMember(int member, int category){
    Category& c = currentCategories[category];
    for(int j = 0; j < c.size; j++){
        if(c.members[j] == member){
            c.size--;
            c.members.erase(c.members.begin() + j);
            if(c.size == 0){
                removeCategory(category);
                return category;
            }
            return -1;
        }
    }
    return -2;
}

int DirichletProcessPrior::popBackCategory(int category){
    Category& c = currentCategories[category];
    c.size--;
    c.members.pop_back();
    if(c.size == 0){
        removeCategory(category);
        return 1;
    }
    return -1;
}

void DirichletProcessPrior::assignMember(int member, int category){
    currentCategories[category].size++;
    currentCategories[category].members.push_back(member);
}

void DirichletProcessPrior::accept() {
    if(moveChoice != -1){
        if(moveChoice == 0){
            betaAcceptCount += 1;
        }
        else if(moveChoice == 1){
            omegaAcceptCount += 1;
        }
        else if(moveChoice == 2){
            rAcceptCount += 1;
        }
    }

    moveChoice = -1;

    for(int i = 0; i < currentCategories.size(); i++){
        currentCategories[i].dirty = false;
    }

    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;
    oldAssignments = assignments;
}

void DirichletProcessPrior::reject() {
    currentCategories = oldCategories;
    currentLnPrior = oldLnPrior;
    assignments = oldAssignments;

    moveChoice = -1;
}

double DirichletProcessPrior::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double randomMove = rng.uniformRv();
    double hastings = 0.0;

    if(randomMove < 0.25) { // Scale Random Beta
        moveChoice = 0;
        betaCount += 1;

        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());

        this->dirty();
        currentCategories[randomCategory].dirty = true;

        double logB = std::log(currentCategories[randomCategory].beta);

        double proposedLogB = logB + betaDelta * Probability::Normal::rv(&rng);
        double proposedB = std::exp(proposedLogB);

        hastings = 0.0;

        currentCategories[randomCategory].beta = proposedB;
    }
    else if(randomMove < 0.50) { // Scale Random Omega
        moveChoice = 1;
        omegaCount += 1;

        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());

        this->dirty();
        currentCategories[randomCategory].dirty = true;

        double logO = std::log(currentCategories[randomCategory].omega);

        double proposedLogO = logO + omegaDelta * Probability::Normal::rv(&rng);
        double proposedO = std::exp(proposedLogO);

        hastings = 0.0;

        currentCategories[randomCategory].omega = proposedO;
    }
    else if(randomMove < 0.75){
        moveChoice = 2;
        rCount += 1;

        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());

        this->dirty();
        currentCategories[randomCategory].dirty = true;

        double logR = std::log(currentCategories[randomCategory].r);

        double proposedLogR = logR + rDelta * Probability::Normal::rv(&rng);
        double proposedR = std::exp(proposedLogR);

        hastings = 0.0;

        currentCategories[randomCategory].r = proposedR;
    }
    else{ // Gibbs sample according to the numGibbsUpdate option
        hastings = INFINITY;
        this->dirty();
        for(int n = 0; n < numGibbsUpdate; n++) {
            int randomMember = (int)(rng.uniformRv() * numMembers);

            int assignment = assignments[randomMember];
            int deleted = unassignMember(randomMember); // This will also delete the group if empty
            if(deleted >= 0){
                model->getTransitionProbability()->deleteQ(deleted);
            }

            std::vector<double> conditionalL;
            int numCats = currentCategories.size();

            tf::Taskflow taskflow;

            for(int i = 0; i < numCats; i++){
                conditionalL.push_back(0.0);

                taskflow.emplace([this, &conditionalL, i, randomMember](){
                    double likelihood = model->testCategory(randomMember, i, false);
                    conditionalL[i] = likelihood + std::log(currentCategories[i].size);
                });
            }

            std::vector<double> omegaVec;
            std::vector<double> betaVec;
            std::vector<double> rVec;
            double alphaSplit = std::log(alpha/10);

            model->getTransitionProbability()->allocateQ(numCats + 10);
            for(int i = 0; i < 10; i++){
                conditionalL.push_back(0.0);
                double newOmega = Probability::Exponential::rv(&rng, omegaLambda);
                double newBeta = Probability::Exponential::rv(&rng, omegaLambda);
                double newR = Probability::Exponential::rv(&rng, rLambda);
                addCategory(newOmega, newBeta, newR);
                omegaVec.push_back(newOmega);
                betaVec.push_back(newBeta);
                rVec.push_back(newR);

                taskflow.emplace([this, &conditionalL, i, randomMember, numCats, alphaSplit](){
                    double likelihood = model->testCategory(randomMember, numCats+i, true);
                    conditionalL[numCats + i] = likelihood + alphaSplit;
                });
            }

            executor.run(taskflow).wait();

            for(int i = 0; i < 10; i++)
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
            for(int i = 0; i < conditionalL.size(); i++){
                total += conditionalL[i];
                if(total > categoryDraw){
                    if(i < numCats) { //It already exists
                        assignMember(randomMember, i);
                        model->getTransitionProbability()->deleteNQ(10);
                        model->regenerateLikelihood(randomMember, i, false);
                    }
                    else {
                        int adjustedIndex = i - numCats;
                        addCategory(omegaVec[adjustedIndex], betaVec[adjustedIndex], rVec[adjustedIndex]);
                        assignMember(randomMember, numCats);
                        model->getTransitionProbability()->deleteNQ(9);
                        model->regenerateLikelihood(randomMember, numCats, true);
                    }
                    break;
                }
            }
        }

        int numCats = currentCategories.size();
        for(int i = 0; i < numCats; i++)
            for(int m : currentCategories[i].members)
                assignments[m] = i;
    }

    regeneratePrior();

    return hastings;
}

void DirichletProcessPrior::tune(){
    double betaRate = (double)betaAcceptCount/(double)betaCount;

    if ( betaRate > 0.44 ) {
        betaDelta *= (1.0 + ((betaRate-0.44)/0.766));
    }
    else {
        betaDelta /= (2.0 - betaRate/0.44);
    }
    betaAcceptCount = 0;
    betaCount = 0;

    double rRate = (double)rAcceptCount/(double)rCount;

    if ( rRate > 0.44 ) {
        rDelta *= (1.0 + ((rRate-0.44)/0.766));
    }
    else {
        rDelta /= (2.0 - rRate/0.44);
    }
    rAcceptCount = 0;
    rCount = 0;


    double omegaRate = (double)omegaAcceptCount/(double)omegaCount;

    if ( omegaRate > 0.44 ) {
        omegaDelta *= (1.0 + ((omegaRate-0.44)/0.766));
    }
    else {
        omegaDelta /= (2.0 - omegaRate/0.44);
    }
    omegaAcceptCount = 0;
    omegaCount = 0;
}