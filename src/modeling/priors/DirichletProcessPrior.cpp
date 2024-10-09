#include "DirichletProcessPrior.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Msg.hpp"
#include <cmath>
#include <iostream>

DirichletProcessPrior::DirichletProcessPrior(int size, double a) : alpha(a), numMembers(size), currentLnPrior(0.0), oldLnPrior(0.0) {
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

    dirty();
    regenerate();
    accept();
}

DirichletProcessPrior::~DirichletProcessPrior(void) {}

int DirichletProcessPrior::getCategorySize(int index) {
    return currentCategories[index].size;
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

void DirichletProcessPrior::setCategoryOmega(int index, double value){
    currentCategories[index].omega = value;
}

void DirichletProcessPrior::setCategoryBeta(int index, double value){
    currentCategories[index].beta = value;
}


double DirichletProcessPrior::getCategoryOmega(int index){
    return currentCategories[index].omega;
}

double DirichletProcessPrior::getCategoryBeta(int index){
    return currentCategories[index].beta;
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

void DirichletProcessPrior::regenerate() {

    if(this->isDirty()){
        for(int i = 0; i < currentCategories.size(); i++)
            for(int m : currentCategories[i].members)
                assignments[m] = i;
        
        currentLnPrior = 0.0;
        for(int& c : assignments){
            currentLnPrior += -2 * std::log(1.0 + currentCategories[c].omega);
            currentLnPrior += Probability::Beta::lnPdf(1, 1, currentCategories[c].beta);
        }
    }
}

void DirichletProcessPrior::sample() {}
