#include "DirichletProcessPrior.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp" // Annoying circular dependency... It is what is is for now...
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/model/ConditionalLikelihood.hpp"
#include "core/Probability.hpp"
#include "core/Msg.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include "core/Math.hpp"
#include <iostream>
#include <algorithm>

DirichletProcessPrior::DirichletProcessPrior(int size, Settings s) : 
                                             alpha(s.dppAlpha), omegaLambda(s.omegaLambda), 
                                             numMembers(size), currentLnPrior(0.0), numGibbsUpdate(s.numGibbsUpdate), 
                                             model(nullptr), omegaDelta(0.25), assignments(size, -1), omegaAcceptCount(0),
                                             omegaCount(0), moveChoice(-1) {}

DirichletProcessPrior::~DirichletProcessPrior() {
    
}

void DirichletProcessPrior::registerModel(Model* m) {
    model = m;
    double expectedCategories = alpha * std::log(1 + numMembers/alpha);
    std::cout << "Initializing Dirichlet Process With E(Categories) = " << expectedCategories << std::endl;
    double flooredCategories = (double)(int)expectedCategories;
    double quantile = 1.0/flooredCategories * numMembers;
    std::cout << "Binning Sites By Heterogeneity With " << flooredCategories << " Bins" << std::endl;

    std::vector<double> heterogeneity;
    std::vector<int> aaMap = {8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 19, 19, 15, 15, 15, 15, 1, 18, 1, 9, 4, 9, 4};  

    ConditionalLikelihood* condL = model->getConditionalLikelihood();
    int numTaxa = model->getNumTaxa();

    for(int i = 0; i < numMembers; i++){
        std::set<int> seenAA;
        for(int j = 0; j < numTaxa; j++){
            double* p = (*condL)(j, 0, 0) + i * 61;
            for(int k = 0; k < 61; k++){
                if(p[k] == 1.0){
                    seenAA.insert(aaMap[k]);
                    break;
                }
            }
        }
        heterogeneity.push_back(seenAA.size());
    }

    std::vector<double> sortedVector(heterogeneity);
    std::sort(sortedVector.begin(), sortedVector.end());

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    for(int i = 0; i < flooredCategories; i++){
        double newOmega1 = Probability::Exponential::rv(&rng, omegaLambda);
        double newOmega2 = Probability::Exponential::rv(&rng, omegaLambda);
        Category newCat = {newOmega1, newOmega2, 0, {}, false}; //Set false so dirty ones aren't copied into "old"
        currentCategories.push_back(newCat);
    }

    std::vector<std::pair<double, int>> indexedHeterogeneity;
    for (int i = 0; i < heterogeneity.size(); i++) {
        indexedHeterogeneity.emplace_back(heterogeneity[i], i);
    }

    std::sort(indexedHeterogeneity.begin(), indexedHeterogeneity.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.first < b.first;
            });

    std::vector<int> sortedIndices;
    for (const auto& pair : indexedHeterogeneity) {
        sortedIndices.push_back(pair.second);
    }

    std::vector<int> thresholds;
    for (int q = 1; q <= flooredCategories; q++) {
        int thresholdIndex = std::min((int)(q * quantile), (int)(sortedIndices.size() - 1));
        thresholds.push_back(thresholdIndex);
    }

    for (int i = 0; i < sortedIndices.size(); i++) {
        for (int c = 0; c < thresholds.size(); c++) {
            if (i < thresholds[c]) {
                assignments[sortedIndices[i]] = c;
                assignMember(sortedIndices[i], c);
                currentCategories[c].members.push_back(sortedIndices[i]);
                currentCategories[c].size++;
                break;
            }
        }
        if (assignments[sortedIndices[i]] == -1) {
            int finalIndex = thresholds.size() - 1;
            assignments[sortedIndices[i]] = finalIndex;
            currentCategories[finalIndex].members.push_back(sortedIndices[i]);
            currentCategories[finalIndex].size++;
        }
    }

    regeneratePrior();

    oldCategories = currentCategories;
    oldLnPrior = currentLnPrior;

    this->dirty();
    for(Category& c : currentCategories)
        c.dirty = true;
}

void DirichletProcessPrior::regeneratePrior(){
    int numCats = currentCategories.size();
    
    currentLnPrior = std::log(alpha) * numCats;

    for(Category& c : currentCategories) {
        currentLnPrior += Math::lnFactorial(c.size - 1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega1);
        currentLnPrior += Probability::Exponential::lnPdf(omegaLambda, c.omega2);
    }
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
    if(moveChoice == 0){
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

double DirichletProcessPrior::updateOmega() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    this->dirty();
    double hastings = 0.0;
    
    moveChoice = 0;
    omegaCount += 1;

    int randomCategory = (int)(rng.uniformRv() * currentCategories.size());
    currentCategories[randomCategory].dirty = true;
    int randomOmega = (int)(rng.uniformRv() * 2);

    if(randomOmega == 0){
        double currentV1 = currentCategories[randomCategory].omega1;
        double scale1 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV1 = currentV1 * scale1;

        currentCategories[randomCategory].omega1 = newV1;
        hastings = std::log(scale1);
    }
    else{
        double currentV2 = currentCategories[randomCategory].omega2;
        double scale2 = std::exp(omegaDelta * (rng.uniformRv() - 0.5));
        double newV2 = currentV2 * scale2;

        currentCategories[randomCategory].omega2 = newV2;
        hastings = std::log(scale2);
    }

    regeneratePrior();

    return hastings;
}

double DirichletProcessPrior::updateDPP(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

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

        std::vector<double> omega1Vec;
        std::vector<double> omega2Vec;
        double alphaSplit = std::log(alpha/5);

        model->getTransitionProbability()->allocateQ(numCats + 5);
        for(int i = 0; i < 5; i++){
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

        for(int i = 0; i < 5; i++)
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
                    model->getTransitionProbability()->deleteNQ(5);
                }
                else {
                    addCategory(omega1Vec[i - numCats], omega2Vec[i - numCats]);
                    assignMember(randomMember, numCats);
                    model->getTransitionProbability()->deleteNQ(4);
                    model->regenerateTransitionProbs(randomMember, numCats);
                }
                break;
            }
        }
    }

    regeneratePrior();

    currentCategories.shrink_to_fit();

    return INFINITY;
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