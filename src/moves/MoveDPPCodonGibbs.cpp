#include "MoveDPPCodonGibbs.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "modeling/likelihoods/MultiMatrixTransitionProbability.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

MoveDPPCodonGibbs::MoveDPPCodonGibbs(CodonMMPhyloCTMC* l, DirichletProcessPrior* d) : likelihood(l), dpp(d), currentMember(0) {}
  
double MoveDPPCodonGibbs::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    int deleted = dpp->unassignMember(currentMember); // This will also delete the group if empty
    if(deleted >= 0){
        likelihood->getTransitionProbability()->deleteQ(deleted);
    }

    dpp->dirty();

    std::vector<double> conditionalL;
    int numCats = dpp->getNumCategories();

    for(int i = 0; i < numCats; i++){
        dpp->assignMember(currentMember, i);
        conditionalL.push_back(likelihood->regenerateAtSite(currentMember, i, false) + std::log(dpp->getCategorySize(i)));
        dpp->popBackCategory(i);
    }

    std::vector<double> newValues;
    double alphaSplit = std::log(dpp->getAlpha()/5);

    for(int i = 0; i < 5; i++){
        double newVal = Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1);
        newValues.push_back(newVal);
        dpp->addCategory(newVal);
        dpp->assignMember(currentMember, numCats);
        conditionalL.push_back(likelihood->regenerateAtSite(i, numCats, true) + alphaSplit);
        dpp->popBackCategory(numCats);
        
    }


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
                dpp->assignMember(currentMember, i);
                likelihood->getTransitionProbability()->deleteQ(numCats); //Remove extra
            }
            else {
                dpp->addCategory(newValues[i - numCats]);
                dpp->assignMember(currentMember, numCats);
            }
            break;
        }
    }

    if(currentMember < dpp->getNumMembers() - 1)
        currentMember++;
    else
        currentMember = 0;

    return INFINITY;
}

void MoveDPPCodonGibbs::tune(){
    acceptedSinceTune = 0;
    countSinceTune = 0;
}