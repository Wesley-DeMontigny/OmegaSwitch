#include "TreeParameter.hpp"
#include "core/RandomVariable.hpp"
#include "TreeObject.hpp"
#include "core/Probability.hpp"
#include "core/Alignment.hpp"
#include "Node.hpp"
#include <cmath>

TreeParameter::TreeParameter(Alignment* aln, std::string newick, double l) : lambda(l), currentPrior(0.0), oldPrior(0.0), 
                                                         branchDelta(0.25), moveChoice(-1), branchCount(0), branchAcceptCount(0), 
                                                         treeCount(0), treeAcceptCount(0), treeDelta(0.25) {
    fixedTree = newick != "";
    if(!fixedTree)
        trees[0] = new TreeObject(aln);
    else
        trees[0] = new TreeObject(newick, aln->getTaxaNames());

    if(!fixedTree){
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
        for(Node* n : nodes) {
            if(n != trees[0]->getRoot()) {
                trees[0]->setBranchLength(n, Probability::Exponential::rv(&rng, lambda));
            }
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
        branchAcceptCount += 1;
    }
    else if(moveChoice == 1){
        treeAcceptCount += 1;
    }

    moveChoice = -1;
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
    
    if(randomMove < 0.9){
        if(!fixedTree){
            moveChoice = 0;
            branchCount += 1;
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

            double scale = std::exp(branchDelta * (rng.uniformRv() - 0.5));

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
        else {
            moveChoice = 0;
            branchCount += 1;
            std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
            Node* root = trees[0]->getRoot();

            Node* p = nullptr;
            do{
                p = nodes[(int)(rng.uniformRv() * nodes.size())];
            }
            while(p == root);

            if(p->getIsTip() == true){
                double currentV = trees[0]->getBranchLength(p);
                double scale = std::exp(branchDelta * (rng.uniformRv() - 0.5));
                double newV = currentV * scale;
                trees[0]->setBranchLength(p, newV);
                p->setNeedsTPUpdate(true);

                Node* q = p;
                do{
                    if(q->getIsTip() == false)
                        q->setNeedsCLUpdate(true);
                    
                    q = q->getAncestor();
                } 
                while(q != root);
                root->setNeedsCLUpdate(true);

                this->dirty();

                hastings = std::log(scale);
            }
            else {
                double currentV = trees[0]->getBranchLength(p);
                double scale = std::exp(branchDelta * (rng.uniformRv() - 0.5));
                double newV = currentV * scale;
                trees[0]->setBranchLength(p, newV);
                p->setNeedsTPUpdate(true);

                for(Node* n : p->getNeighbors()){
                    if(n != p->getAncestor()){
                        double currentNLength = trees[0]->getBranchLength(n);
                        double newLength = currentNLength * scale;
                        trees[0]->setBranchLength(n, newLength);
                        n->setNeedsTPUpdate(true);
                        if(n->getIsTip() == false)
                            n->setNeedsCLUpdate(true);
                    }
                }

                Node* q = p;
                do{
                    if(q->getIsTip() == false)
                        q->setNeedsCLUpdate(true);
                    
                    q = q->getAncestor();
                } 
                while(q != root);
                root->setNeedsCLUpdate(true);

                this->dirty();

                hastings = 3 * std::log(scale);
            }
        }
    }
    else {
        moveChoice = 1;
        treeCount += 1;
        std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
        trees[0]->updateAll();
        this->dirty();
        
        Node* root = trees[0]->getRoot();
        double scale = std::exp(treeDelta * (rng.uniformRv() - 0.5));

        for(Node* n : trees[0]->getPostOrderSeq()){
            if(n != root){
                trees[0]->setBranchLength(n, trees[0]->getBranchLength(n) * scale);
            }
        }

        hastings = (trees[0]->getNumNodes()-1) * std::log(scale);
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
    double rate1 = (double)branchAcceptCount/(double)branchCount;

    if ( rate1 > 0.25 ) {
        branchDelta *= (1.0 + ((rate1-0.25)/0.766));
    }
    else {
        branchDelta /= (2.0 - rate1/0.25);
    }
    branchAcceptCount = 0;
    branchCount = 0;

    double rate2 = (double)treeAcceptCount/(double)treeCount;

    if ( rate2 > 0.25 ) {
        treeDelta *= (1.0 + ((rate2-0.25)/0.766));
    }
    else {
        treeDelta /= (2.0 - rate2/0.25);
    }
    treeAcceptCount = 0;
    treeCount = 0;
}

double TreeParameter::lnPrior() {
    return currentPrior;
}