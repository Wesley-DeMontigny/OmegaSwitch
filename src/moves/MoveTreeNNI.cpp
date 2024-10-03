#include "MoveTreeNNI.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "MoveScheduler.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/parameters/trees/Node.hpp"

MoveTreeNNI::MoveTreeNNI(TreeParameter* t) : param(t) {}
        
double MoveTreeNNI::update(){

    TreeObject* tree = param->getTree();
    std::vector<Node*> nodes = tree->getPostOrderSeq();
    Node* root = tree->getRoot();

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Node* p = nullptr;
    do{
        p = nodes[(int)(rng.uniformRv() * nodes.size())];
    }
    while(p == root || p->getIsTip() == true);

    Node* a = p->getAncestor();

    std::set<Node*> neighbors1 = p->getNeighbors();
    neighbors1.erase(a);//Exclude a
    Node* n1 = Node::chooseNodeFromSet(neighbors1);

    std::set<Node*> neighbors2 = a->getNeighbors();
    neighbors2.erase(p);//Don't select p
    Node* n2 = Node::chooseNodeFromSet(neighbors2);

    n1->addNeighbor(a);
    a->addNeighbor(n1);
    n1->setAncestor(a);
    
    if(n2 != a->getAncestor()){
        n2->addNeighbor(p);
        p->addNeighbor(n2);
        n2->setAncestor(p);
    }
    else{//If this isn't the case then we need to swap around the tree
        n2->addNeighbor(p);
        p->addNeighbor(n2);
        p->setAncestor(n2);
        a->setAncestor(p);
    }

    //Remove old connections
    n1->removeNeighbor(p);
    p->removeNeighbor(n1);
    n2->removeNeighbor(a);
    a->removeNeighbor(n2);

    Node* q = nullptr;
    if(n2 != a->getAncestor())
        q = p;
    else{
        q = a;
        a->setNeedsTPUpdate(true);//Because TP is defined as the branch going to the descendent...
        p->setNeedsTPUpdate(true);//Flipping this internal branch causes a change in the TP!
    }

    do{
        q->setNeedsCLUpdate(true);
        q = q->getAncestor();
    }
    while(q != root);
    root->setNeedsCLUpdate(true);

    tree->initPostOrder();

    param->dirty();

    return 0.0;
}

void MoveTreeNNI::tune(){
    acceptedSinceTune = 0;
    countSinceTune = 0;
}