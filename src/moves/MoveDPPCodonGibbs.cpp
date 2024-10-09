#include "MoveDPPCodonGibbs.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "modeling/likelihoods/TransitionProbability.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

MoveDPPCodonGibbs::MoveDPPCodonGibbs(PhyloCTMC* l, DirichletProcessPrior* d) : likelihood(l), dpp(d) {}
  
double MoveDPPCodonGibbs::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    int currentMember = (int)(rng.uniformRv() * dpp->getNumMembers());

    int assignment = dpp->getAssinments()[currentMember];
    int deleted = dpp->unassignMember(currentMember); // This will also delete the group if empty
    if(deleted >= 0){
        likelihood->getTransitionProbability()->deleteQ(deleted);
    }

    std::vector<double> conditionalL;
    int numCats = dpp->getNumCategories();

    for(int i = 0; i < numCats; i++){
        dpp->assignMember(currentMember, i);
        conditionalL.push_back(likelihood->regenerateIntoSiteBuffer(currentMember, i, false) + std::log(dpp->getCategorySize(i)));
        dpp->popBackCategory(i);
    }

    std::vector<double> omegaVec;
    std::vector<double> betaVec;
    double alphaSplit = std::log(dpp->getAlpha()/5);

    for(int i = 0; i < 5; i++){
        double newOmega = Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1);
        double newBeta = Probability::Beta::rv(&rng, 1, 1);
        omegaVec.push_back(newOmega);
        betaVec.push_back(newBeta);
        dpp->addCategory(newOmega, newBeta);
        dpp->assignMember(currentMember, numCats);
        conditionalL.push_back(likelihood->regenerateIntoSiteBuffer(i, numCats, true) + alphaSplit);
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
    int newAssignment = -1;

    total = 0.0;
    for(int i = 0; i < conditionalL.size(); i++){
        total += conditionalL[i];
        if(total > categoryDraw){
            if(i < numCats) { //It already exists
                dpp->assignMember(currentMember, i);
                newAssignment = i;
                likelihood->getTransitionProbability()->deleteQ(numCats); //Remove extra
            }
            else {
                dpp->addCategory(omegaVec[i - numCats], betaVec[i - numCats]);
                newAssignment = numCats;
                dpp->assignMember(currentMember, numCats);
            }
            break;
        }
    }

    if(newAssignment != assignment){
        dpp->dirty();
    }

    return INFINITY;
}

void MoveDPPCodonGibbs::tune(){
    acceptedSinceTune = 0;
    countSinceTune = 0;
}