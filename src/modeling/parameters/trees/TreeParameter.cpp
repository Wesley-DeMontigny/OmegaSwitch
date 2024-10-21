#include "TreeParameter.hpp"
#include "core/RandomVariable.hpp"
#include "TreeObject.hpp"
#include "core/Probability.hpp"
#include "Node.hpp"

TreeParameter::TreeParameter(Alignment* aln, double l) : lambda(l), currentPrior(0.0), oldPrior(0.0), 
                                                         localDelta(std::log(4.0)), moveChoice(-1), 
                                                         localCount(0), localAcceptCount(0) {
    trees[0] = new TreeObject(aln);

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
    for(Node* n : nodes) {
        if(n != trees[0]->getRoot()) {
            trees[0]->setBranchLength(n, Probability::Exponential::rv(&rng, lambda));
        }
    }

    trees[1] = new TreeObject(*trees[0]);

    std::vector<double> values = trees[0]->getBranchLengths();
    double totalLength = 0.0;
    for(double val : values){
        totalLength += val;
    }
    currentPrior = Probability::Gamma::lnPdf(values.size(), lambda, totalLength);
    oldPrior = currentPrior;

    dirty();
}

TreeParameter::~TreeParameter(){
    delete trees[0];
    delete trees[1]; 
}

void TreeParameter::accept(){
    *trees[1] = *trees[0];
    oldPrior = currentPrior;

    if(moveChoice == 0){
        localAcceptCount += 1;
        moveChoice = -1;
    }
}

void TreeParameter::reject(){
    *trees[0] = *trees[1];
    currentPrior = oldPrior;

    moveChoice = -1;
}


double TreeParameter::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double randomMove = rng.uniformRv();

    double hastings = 0.0;

    //Local Move
    if(randomMove < 0.33){
        moveChoice = 0;
        localCount += 1;
        TreeObject* tree = trees[0];
        std::vector<Node*> nodes = tree->getPostOrderSeq();
        Node* root = tree->getRoot();

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

        double scale = std::exp(localDelta * (rng.uniformRv() - 0.5));

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
        this->dirty();

        hastings = 3 * std::log(scale);
    }
    else { //NNI
        TreeObject* tree = trees[0];
        std::vector<Node*> nodes = tree->getPostOrderSeq();
        Node* root = tree->getRoot();

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
        this->dirty();

        hastings = 0.0;
    }

    std::vector<double> values = trees[0]->getBranchLengths();
    double totalLength = 0.0;
    for(double val : values){
        totalLength += val;
    }
    currentPrior = Probability::Gamma::lnPdf(values.size(), lambda, totalLength);

    return hastings;
}

void TreeParameter::tune() {
    double rate = (double)localAcceptCount/(double)localCount;

    if ( rate > 0.44 ) {
        localDelta *= (1.0 + ((rate-0.44)/0.766));
    }
    else {
        localDelta /= (2.0 - rate/0.44);
    }
    localAcceptCount = 0;
    localCount = 0;
}

double TreeParameter::lnPrior() {
    return currentPrior;
}