#include "MoveTreeLocal.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "MoveScheduler.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <iostream>
#include <cmath>

MoveTreeLocal::MoveTreeLocal(TreeParameter* t) : param(t) {}
        
/*
See Larget and Simon (1999) for original description of the move and Mark Holder
et al. (2005) for the correction on the Hastings ratio.

The move is basically as follows:
1. Pick an internal branch and two neighbor nodes of each of the nodes of the branch
2. Scale the length of that whole path like you would normally do.
3. Randomly draw a position on that length to pick off the subtree originally connected
   to one of the internal nodes and graft it there.
*/
double MoveTreeLocal::update(){

    TreeObject* tree = param->getTree();
    std::vector<Node*> nodes = tree->getPostOrderSeq();
    Node* root = tree->getRoot();

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Node* u = nullptr;
    do{
        u = nodes[(int)(rng.uniformRv() * nodes.size())];
    }
    while(u == root || u->getIsTip() == true);
    Node* v = u->getAncestor();

    std::set<Node*> neighbors1 = u->getNeighbors();
    neighbors1.erase(v);//Exclude v
    Node* a = Node::chooseNodeFromSet(neighbors1);

    std::set<Node*> neighbors2 = v->getNeighbors();
    neighbors2.erase(u);//Don't select u
    Node* c = Node::chooseNodeFromSet(neighbors2);

    double scale = std::exp(delta * (rng.uniformRv() - 0.5));

    double paths[3];
    paths[0] = tree->getBranchLength(u) * scale;
    paths[1] = tree->getBranchLength(a) * scale;
    Node* b3 = nullptr;
    if(c != v->getAncestor())
        b3 = c;
    else
        b3 = v;
    paths[2] = tree->getBranchLength(b3) * scale;

    double totalPath = paths[0] + paths[1] + paths[2];

    std::vector<Node*> nodeSet = {a, b3};
    int pick = (int)(rng.uniformRv() * 2);
    double randomLoc = rng.uniformRv() * totalPath;

    //The pick decides the oritentation of the path
    if(randomLoc <= totalPath - paths[pick]){
        tree->setBranchLength(nodeSet[pick], randomLoc);
        tree->setBranchLength(u, totalPath - paths[pick] - randomLoc);
        tree->setBranchLength(nodeSet[pick ^ 1], paths[pick]);
    }
    else{
        u->removeNeighbor(a);
        a->removeNeighbor(u);
        v->removeNeighbor(c);
        c->removeNeighbor(v);

        v->addNeighbor(a);
        a->addNeighbor(v);
        u->addNeighbor(c);
        c->addNeighbor(u);
        tree->setBranchLength(nodeSet[pick ^ 1], totalPath - randomLoc);
        tree->setBranchLength(u, randomLoc - (totalPath - paths[pick]));
        tree->setBranchLength(nodeSet[pick], totalPath - paths[pick]);

        //Rooting logic
        if(v->getAncestor() == c){
            u->setAncestor(c);
            v->setAncestor(u);
            a->setAncestor(v);
        }
        else{
            c->setAncestor(u);
            u->setAncestor(v);
            a->setAncestor(v);
        }
    }


    u->setNeedsTPUpdate(true);
    v->setNeedsTPUpdate(true);
    a->setNeedsTPUpdate(true);
    c->setNeedsTPUpdate(true);

    //The updating gets a little awkward because we don't really know the branching here.
    Node* q = v;
    if(v->getAncestor() != u)
        q = u;

    do{
        if(q->getIsTip() == false)
            q->setNeedsCLUpdate(true);
        q = q->getAncestor();
    }
    while(q != root);
    root->setNeedsCLUpdate(true);

    tree->initPostOrder();

    param->dirty();

    return std::pow(scale, 3);
}

void MoveTreeLocal::tune(){
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