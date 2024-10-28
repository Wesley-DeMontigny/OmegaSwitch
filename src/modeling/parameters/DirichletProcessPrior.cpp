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
                                             oldLnPrior(0.0), numGibbsUpdate(numGibbs), model(nullptr),
                                             moveChoice(-1), omegaCount(0), omegaAcceptCount(0), omegaDelta(std::log(2)),
                                             betaCount(0), betaAcceptCount(0), betaAlpha(1.0) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //The expected category number can be roughly computed as n = a * ln(1 + c/a)
    for(int i = 0; i < size; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            Category newCat = {Probability::Exponential::rv(&rng, omegaLambda),
                               Probability::Beta::rv(&rng, 1, 1),
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
        currentLnPrior += Probability::Beta::lnPdf(1, 1, currentCategories[i].beta);
    }

    currentLnPrior -= denominator;
}

void DirichletProcessPrior::removeCategory(int index){
    if(currentCategories[index].size != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    currentCategories.erase(currentCategories.begin() + index);
}

void DirichletProcessPrior::addCategory(double value1, double value2){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Category newCat = {value1, value2, 0, {}, true};
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

    if(randomMove < 0.33) { // Beta Simplex on Random Beta
        moveChoice = 0;
        betaCount += 1;

        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());

        this->dirty();
        currentCategories[randomCategory].dirty = true;

        double betaVal = currentCategories[randomCategory].beta;

        double a = betaAlpha + 1.0;
        double b = (betaAlpha / betaVal) - a + 2.0;
        double newVal = Probability::Beta::rv(&rng, a, b);

        currentCategories[randomCategory].beta = newVal;

        double scalingFactor = (1.0 - newVal)/(1.0 - betaVal);

        double forward = Probability::Beta::lnPdf(a, b, newVal);
        double newA = betaAlpha + 1.0;
        double newB = (betaAlpha / newVal) - a + 2.0;
        double backward = Probability::Beta::lnPdf(newA, newB, betaVal);
        
        hastings = backward - forward;
    }
    else if(randomMove < 0.66) { // Scale Random Omega
        moveChoice = 1;
        omegaCount += 1;

        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());

        this->dirty();
        currentCategories[randomCategory].dirty = true;

        double scale = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        currentCategories[randomCategory].omega *= scale;

        hastings = std::log(scale);
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
            for(int i = 0; i < numCats; i++){
                model->regenerateLikelihood(randomMember, i, false);
                conditionalL.push_back(model->lnLikelihood() + std::log(currentCategories[i].size));
            }

            std::vector<double> omegaVec;
            std::vector<double> betaVec;
            double alphaSplit = std::log(alpha/5);

            addCategory(0, 0);
            for(int i = 0; i < 5; i++){
                double newOmega = Probability::Exponential::rv(&rng, omegaLambda);
                double newBeta = Probability::Beta::rv(&rng, 1, 1);
                currentCategories[numCats].omega = newOmega;
                currentCategories[numCats].beta = newBeta;
                omegaVec.push_back(newOmega);
                betaVec.push_back(newBeta);
                model->regenerateLikelihood(randomMember, numCats, true);
                conditionalL.push_back(model->lnLikelihood() + alphaSplit);
            }
            removeCategory(numCats);


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
                        model->getTransitionProbability()->deleteQ(numCats); //Remove extra
                        model->regenerateLikelihood(randomMember, i, false);
                    }
                    else {
                        addCategory(omegaVec[i - numCats], betaVec[i - numCats]);
                        assignMember(randomMember, numCats);
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
        betaAlpha *= (1.0 + ((betaRate-0.44)/0.766));
    }
    else {
        betaAlpha /= (2.0 - betaRate/0.44);
    }
    betaAcceptCount = 0;
    betaCount = 0;


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