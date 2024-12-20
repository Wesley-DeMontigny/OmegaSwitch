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

DirichletProcessPrior::DirichletProcessPrior(int size, double a, double oL, int numGibbs) : 
                                             alpha(a), omegaLambda(oL), numMembers(size), currentLnPrior(0.0), 
                                             numGibbsUpdate(numGibbs), model(nullptr), omegaDelta(0.15) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    for(int i = 0; i < size; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            Category newCat = {Probability::Exponential::rv(&rng, omegaLambda),
                               Probability::Exponential::rv(&rng, omegaLambda),
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

    //This is just a constant - no need to calculate it
    //denominator = Math::lnGamma(numMembers - alpha);

    int numCats = currentCategories.size();
    for(int i = 0; i < numCats; i++)
        for(int m : currentCategories[i].members)
            assignments[m] = i;
    
    this->dirty();

    regeneratePrior();

    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;
}

DirichletProcessPrior::~DirichletProcessPrior() {
    
}

void DirichletProcessPrior::regeneratePrior(){
    int numCats = currentCategories.size();
    
    currentLnPrior = std::log(alpha) * numCats;

    for(Category& c : currentCategories) {
        currentLnPrior += Math::lnFactorial(c.size - 1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega2);
    }

    //currentLnPrior -= denominator;
}

void DirichletProcessPrior::removeCategory(int index){
    if(currentCategories[index].size != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    currentCategories.erase(currentCategories.begin() + index);
}

void DirichletProcessPrior::addCategory(double omega1, double omega2){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Category newCat = {omega1, omega2, 0, {}};
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
    if(moveChoice != 0){
        omegaAcceptCount += 1;
    }

    moveChoice = -1;

    for(int i = 0; i < currentCategories.size(); i++){
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

double DirichletProcessPrior::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    this->dirty();
    double hastings = 0.0;

    double randomMove = rng.uniformRv();

    if(randomMove < 0.5) { // Scale Random Omega
        moveChoice = 0;
        omegaCount += 1;

        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());
        int randomOmega = (int)(rng.uniformRv() * 2);

        currentCategories[randomCategory].dirty = true;

        if(randomOmega == 0){
            double logO = std::log(currentCategories[randomCategory].omega1);

            double proposedLogO = logO + omegaDelta * Probability::Normal::rv(&rng);
            double proposedO = std::exp(proposedLogO);

            currentCategories[randomCategory].omega1 = proposedO;
        }
        else{
            double logO = std::log(currentCategories[randomCategory].omega2);

            double proposedLogO = logO + omegaDelta * Probability::Normal::rv(&rng);
            double proposedO = std::exp(proposedLogO);

            currentCategories[randomCategory].omega2 = proposedO;
        }
    }
    else{
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

            std::vector<double> omega1Vec;
            std::vector<double> omega2Vec;
            double alphaSplit = std::log(alpha/20);

            model->getTransitionProbability()->allocateQ(numCats + 20);
            for(int i = 0; i < 20; i++){
                conditionalL.push_back(0.0);
                double newOmega1 = Probability::Exponential::rv(&rng, omegaLambda);
                double newOmega2 = Probability::Exponential::rv(&rng, omegaLambda);
                addCategory(newOmega1, newOmega2);
                omega1Vec.push_back(newOmega1);
                omega2Vec.push_back(newOmega2);

                taskflow.emplace([this, &conditionalL, i, randomMember, numCats, alphaSplit](){
                    double likelihood = model->testCategory(randomMember, numCats+i, true);
                    conditionalL[numCats + i] = likelihood + alphaSplit;
                });
            }

            executor.run(taskflow).wait();

            for(int i = 0; i < 20; i++)
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
                        model->getTransitionProbability()->deleteNQ(20);
                        model->regenerateLikelihood(randomMember, i, false);
                    }
                    else {
                        addCategory(omega1Vec[i - numCats], omega2Vec[i - numCats]);
                        assignMember(randomMember, numCats);
                        model->getTransitionProbability()->deleteNQ(19);
                        model->regenerateLikelihood(randomMember, numCats, true);
                    }
                    break;
                }
            }
        }
        hastings = INFINITY;
    }

    regeneratePrior();

    return hastings;
}

void DirichletProcessPrior::tune() {
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