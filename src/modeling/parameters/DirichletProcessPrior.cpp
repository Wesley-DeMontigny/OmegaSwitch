#include "DirichletProcessPrior.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp" // Annoying circular dependency... It is what is is for now...
#include "modeling/model/TransitionProbability.hpp"
#include "core/Probability.hpp"
#include "core/Msg.hpp"
#include <cmath>
#include "core/Math.hpp"
#include <iostream>

DirichletProcessPrior::DirichletProcessPrior(int size, double a, int numGibbs) : alpha(a), numMembers(size), currentLnPrior(0.0), 
                                                                                 oldLnPrior(0.0), numGibbsUpdate(numGibbs), model(nullptr) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //The expected category number can be roughly computed as n = a * ln(1 + c/a)
    for(int i = 0; i < size; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            Category newCat = {Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1),
                               Probability::Beta::rv(&rng, 1, 1),
                               1, {i}};
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

    denominator = 0.0;
    double cp = alpha - 1;
    for(int i = 1; i <= numMembers; i++){
        denominator += std::log(cp + i);
    }


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
    currentLnPrior += Math::lnStirlingFirst(numMembers, numCats);

    for (int i = 0; i < numCats; ++i) {
        currentLnPrior += Math::lnFactorial(currentCategories[i].size - 1);
        
        currentLnPrior += -2 * std::log(1.0 + currentCategories[i].omega);
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

    Category newCat = {value1, value2, 0, {}};
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
    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;
    oldAssignments = assignments;
}

void DirichletProcessPrior::reject() {
    currentCategories = oldCategories;
    currentLnPrior = oldLnPrior;
    assignments = oldAssignments;
}

double DirichletProcessPrior::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double randomMove = rng.uniformRv();
    double hastings = 0.0;

    if(randomMove < 0.33) { // Beta Simplex on Random Beta
        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());
        double betaVal = currentCategories[randomCategory].beta;

        double a = alpha + 1.0;
        double b = (alpha / betaVal) - a + 2.0;
        double newVal = Probability::Beta::rv(&rng, a, b);

        currentCategories[randomCategory].beta = newVal;

        double scalingFactor = (1.0 - newVal)/(1.0 - betaVal);

        double forward = Probability::Beta::lnPdf(a, b, newVal);
        double newA = alpha + 1.0;
        double newB = (alpha / newVal) - a + 2.0;
        double backward = Probability::Beta::lnPdf(newA, newB, betaVal);
        
        hastings = backward - forward;

        regeneratePrior();

        this->dirty();
    }
    else if(randomMove < 0.66) { // Scale Random Omega
        double delta = std::log(4);
        int randomCategory = (int)(rng.uniformRv() * currentCategories.size());

        double scale = std::exp(delta * (rng.uniformRv() - 0.5));
        currentCategories[randomCategory].omega *= scale;

        hastings = std::log(scale);

        regeneratePrior();

        this->dirty();
    }
    else { // Gibbs sample according to the numGibbsUpdate option
        hastings = INFINITY;

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
                assignMember(randomMember, i);
                conditionalL.push_back(model->regenerateIntoLikelihoodBuffer(randomMember, i, false) + std::log(currentCategories[i].size));
                popBackCategory(i);
            }

            std::vector<double> omegaVec;
            std::vector<double> betaVec;
            double alphaSplit = std::log(alpha/5);

            for(int i = 0; i < 5; i++){
                double newOmega = Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1);
                double newBeta = Probability::Beta::rv(&rng, 1, 1);
                omegaVec.push_back(newOmega);
                betaVec.push_back(newBeta);
                addCategory(newOmega, newBeta);
                assignMember(randomMember, numCats);
                conditionalL.push_back(model->regenerateIntoLikelihoodBuffer(i, numCats, true) + alphaSplit);
                popBackCategory(numCats);
                
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
                        assignMember(randomMember, i);
                        model->getTransitionProbability()->deleteQ(numCats); //Remove extra
                        if(i != assignment){
                            model->forceRegenerate(randomMember, i, false);
                        }
                    }
                    else {
                        addCategory(omegaVec[i - numCats], betaVec[i - numCats]);
                        assignMember(randomMember, numCats);
                        model->forceRegenerate(randomMember, i, true);
                    }
                    break;
                }
            }
        }


        int numCats = currentCategories.size();
        for(int i = 0; i < numCats; i++)
            for(int m : currentCategories[i].members)
                assignments[m] = i;
        
        regeneratePrior();
    }

    return hastings;
}