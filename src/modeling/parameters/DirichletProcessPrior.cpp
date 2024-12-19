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
                                             numGibbsUpdate(numGibbs), model(nullptr) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    for(int i = 0; i < size * 2; i++){
        double randomVal = rng.uniformRv();
        double total = alpha/(i + alpha);

        // If new category
        if(total > randomVal){
            Category newCat = {Probability::Exponential::rv(&rng, omegaLambda),
                               1, {i}};
            categories.push_back(newCat);
            continue;
        }

        for(Category &c : categories){
            total += c.size/(i+alpha);

            //If old category
            if(total > randomVal){
                c.size++;
                c.members.push_back(i);
                break;
            }
        }
    }

    for(int i = 0; i < numMembers * 2; i++)
        assignments.push_back(-1);

    //This is just a constant - no need to calculate it
    //denominator = Math::lnGamma(numMembers - alpha);

    for(Category& c : categories)
        for(int m : c.members)
            assignments[m] = c.value;
    
    this->dirty();

    regeneratePrior();
}

DirichletProcessPrior::~DirichletProcessPrior() {
    
}

void DirichletProcessPrior::regeneratePrior(){
    int numCats = categories.size();
    
    currentLnPrior = std::log(alpha) * numCats;

    for(Category& c : categories) {
        currentLnPrior += Math::lnFactorial(c.size - 1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.value);
    }

    //currentLnPrior -= denominator;
}

void DirichletProcessPrior::removeCategory(int index){
    if(categories[index].size != 0)
        Msg::error("Attempt to remove DPP category that still had members!");

    categories.erase(categories.begin() + index);
}

void DirichletProcessPrior::addCategory(double omega){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Category newCat = {omega, 0, {}};
    categories.push_back(newCat);
}

int DirichletProcessPrior::unassignMember(int member){
    for(int i = 0; i < categories.size(); i++){
        Category& c = categories[i];
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
    categories[category].size++;
    categories[category].members.push_back(member);
}

void DirichletProcessPrior::accept() {}

void DirichletProcessPrior::reject() {}

double DirichletProcessPrior::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    for(int n = 0; n < numGibbsUpdate; n++) {
        int randomSite = (int)(rng.uniformRv() * numMembers);
        int randomOmega = (int)(rng.uniformRv() * 2);
        int randomMember = randomSite + randomOmega;

        int assignment = assignments[randomMember];
        int deleted = unassignMember(randomMember); // This will also delete the group if empty

        std::vector<double> conditionalL;
        int numCats = categories.size();

        tf::Taskflow taskflow;

        for(int i = 0; i < numCats; i++){
            conditionalL.push_back(0.0);

            double catOmega = categories[i].value;

            taskflow.emplace([this, &conditionalL, i, randomSite, randomOmega, catOmega](){
                double likelihood = model->testCategory(randomSite, randomOmega, catOmega);
                conditionalL[i] = likelihood + std::log(categories[i].size);
            });
        }

        std::vector<double> omegaVec;
        double alphaSplit = std::log(alpha/5);

        for(int i = 0; i < 5; i++){
            conditionalL.push_back(0.0);
            double newOmega = Probability::Exponential::rv(&rng, omegaLambda);

            omegaVec.push_back(newOmega);

            taskflow.emplace([this, &conditionalL, i, randomSite, randomOmega, newOmega, numCats, alphaSplit](){
                double likelihood = model->testCategory(randomSite, randomOmega, newOmega);
                conditionalL[numCats + i] = likelihood + alphaSplit;
            });
        }

        executor.run(taskflow).wait();

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
                    assignments[randomMember] = categories[i].value;
                    model->regenerateLikelihood(randomSite);
                }
                else {
                    int adjustedIndex = i - numCats;
                    addCategory(omegaVec[adjustedIndex]);
                    assignments[randomMember] = omegaVec[adjustedIndex];
                    assignMember(randomMember, numCats);
                    model->regenerateLikelihood(randomMember);
                }
                break;
            }
        }
    }

    regeneratePrior();

    return INFINITY;
}

void DirichletProcessPrior::tune() {}