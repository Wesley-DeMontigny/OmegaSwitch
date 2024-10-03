#include "MoveDPPCodonGibbs.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "modeling/likelihoods/TransitionProbability.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

MoveDPPCodonGibbs::MoveDPPCodonGibbs(PhyloCTMC* l, DirichletProcessPrior* d) : likelihood(l), dpp(d), currentMember(0) {}
  
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

    std::vector<double> omega1Vec;
    std::vector<double> omega2Vec;
    double alphaSplit = std::log(dpp->getAlpha()/5);

    for(int i = 0; i < 5; i++){
        double newOmega1 = Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1);
        double newOmega2 = Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1);
        omega1Vec.push_back(newOmega1);
        omega2Vec.push_back(newOmega2);
        dpp->addCategory(newOmega1, newOmega2);
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
                dpp->addCategory(omega1Vec[i - numCats], omega2Vec[i - numCats]);
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