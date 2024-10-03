#include "TreePrior.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "core/Probability.hpp"
#include "core/RandomVariable.hpp"
#include <cmath>
#include <iostream>

TreePrior::TreePrior(TreeParameter* t) : currentLnPrior(0.0), oldLnPrior(0.0), lambda(nullptr), tree(t) {
    numBranches = tree->getBranchLengths().size();
    this->dirty();
}

void TreePrior::setExponentialBranchPrior(BasicParameter<double>* l){
    lambda = l;
}

void TreePrior::addMonophyleticConstraint(std::set<int> taxa, double strength){
    //Implement later
}

void TreePrior::accept() {
    oldLnPrior = currentLnPrior;
    tree->accept();
    tree->clean();
    if(lambda != nullptr){
        lambda->accept();
        lambda->clean();
    }
}

void TreePrior::reject() {
    currentLnPrior = oldLnPrior;
    tree->reject();
    tree->clean();
    if(lambda != nullptr){
        lambda->reject();
        lambda->clean();
    }
}

void TreePrior::regenerate(){
    tree->regenerate();

    //This isn't pretty. The tree prior should eventually be split between a branch length and topology prior.
    bool hasLambda = lambda != nullptr;
    if(hasLambda){
        lambda->regenerate();
        if(tree->isDirty() || lambda->isDirty())
            this->dirty();
    }
    else{
        if(tree->isDirty())
            this->dirty();
    }
    

    if(this->isDirty()){
        if(hasLambda){
            std::vector<double> values = tree->getBranchLengths();
            double lambdaVal = lambda->getValue();
            double totalLength = 0.0;
            for(double val : values){
                totalLength += val;
            }
            currentLnPrior = Probability::Gamma::lnPdf(numBranches, lambdaVal, totalLength);
        }
    }
}

void TreePrior::sample(){
    if(lambda != nullptr){
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        TreeObject* t = tree->getTree();
        double lambdaVal = lambda->getValue();
        std::vector<Node*> nodes = t->getPostOrderSeq();
        for(Node* n : nodes){
            if(n != t->getRoot()){
                t->setBranchLength(n, Probability::Exponential::rv(&rng, lambdaVal));
            }
        }
    }
    tree->accept();
    tree->regenerate();
}