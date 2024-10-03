#include "PosteriorNode.hpp"
#include <iostream>

PosteriorNode::PosteriorNode(LikelihoodNode* lN, std::vector<PriorNode*> pN) : likelihood(lN), priors(pN), currentLnPosterior(0.0), oldLnPosterior(0.0) {
    this->dirty();
}

void PosteriorNode::accept() {
    oldLnPosterior = currentLnPosterior;

    likelihood->accept();
    likelihood->clean();
    for(PriorNode* p : priors){
        p->accept();
        p->clean();
    }
}

void PosteriorNode::reject() {
    currentLnPosterior = oldLnPosterior;

    likelihood->reject();
    likelihood->clean();
    for(PriorNode* p : priors){
        p->reject();
        p->clean();
    }
}

void PosteriorNode::regenerate(){
    likelihood->regenerate();
    bool pDirty = false;
    for(PriorNode* p : priors){
        p->regenerate();
        if(p->isDirty())
            pDirty = true;
    }
    if(pDirty || likelihood->isDirty())
        this->dirty();

    if(this->isDirty()){
        double pTotal = 0.0;
        for(PriorNode* p : priors)
            pTotal += p->lnPrior();
        currentLnPosterior = likelihood->lnLikelihood() + pTotal;
    }
}