#include "MoveScaleBranch.hpp"
#include "MoveScheduler.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <iostream>

MoveScaleBranch::MoveScaleBranch(TreeParameter* t) : param(t), delta(std::log(4.0)) {}
        
double MoveScaleBranch::update(){

    TreeObject* tree = param->getTree();
    std::vector<Node*> nodes = tree->getPostOrderSeq();
    Node* root = tree->getRoot();

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Node* p = nullptr;
    do{
        p = nodes[(int)(rng.uniformRv() * nodes.size())];
    }
    while(p == root);

    double currentV = tree->getBranchLength(p);
    double scale = std::exp(delta * (rng.uniformRv() - 0.5));
    double newV = currentV * scale;
    tree->setBranchLength(p, newV);
    p->setNeedsTPUpdate(true);

    Node* q = p;
    do{
        if(q->getIsTip() == false)//The conditional likelihoods at the tips should never change
            q->setNeedsCLUpdate(true);
        
        q = q->getAncestor();
    } 
    while(q != root);
    root->setNeedsCLUpdate(true);

    param->dirty();

    return scale;
}

void MoveScaleBranch::tune(){
    double rate = (double)acceptedSinceTune/(double)countSinceTune;

    if ( rate > 0.44 ) {
        delta *= (1.0 + ((rate-0.44)/0.766));
    }
    else {
        delta /= (2.0 - rate/0.44);
    }

    acceptedSinceTune = 0;
    countSinceTune = 0;
}