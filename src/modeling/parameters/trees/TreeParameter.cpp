#include "core/Alignment.hpp"
#include "core/Probability.hpp"
#include "core/RandomVariable.hpp"
#include "Node.hpp"
#include "TreeObject.hpp"
#include "TreeParameter.hpp"
#include <cmath>

TreeParameter::TreeParameter(Alignment& aln, std::string& newick, double l) : lambda(l), currentPrior(0.0), oldPrior(0.0), 
                                                         branchDelta(0.5), moveChoice(TreeMoves::NO_MOVE), branchCount(0), 
                                                         branchAcceptCount(0), treeCount(0), treeAcceptCount(0), treeAlpha(10000) {
    fixedTree = newick != "";
    if(!fixedTree)
        trees[0] = new TreeObject(aln, false);
    else
        trees[0] = new TreeObject(newick, aln.getTaxaNames());

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

TreeParameter::TreeParameter(TreeObject& tree, double lambda) : lambda(lambda), currentPrior(0.0), oldPrior(0.0), branchDelta(0.5), 
                                                                moveChoice(TreeMoves::NO_MOVE), branchCount(0), branchAcceptCount(0), 
                                                                treeCount(0), treeAcceptCount(0), treeAlpha(10000), fixedTree(true) {

    trees[0] = new TreeObject(tree);

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

    if(moveChoice == TreeMoves::BRANCH_MOVE){
        branchAcceptCount += 1;
    }
    else if(moveChoice == TreeMoves::TREE_MOVE){
        treeAcceptCount += 1;
    }

    moveChoice = TreeMoves::NO_MOVE;
}

void TreeParameter::reject(){
    *trees[0] = *trees[1];
    currentPrior = oldPrior;

    moveChoice = TreeMoves::NO_MOVE;
}


double TreeParameter::update() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double randomMove = rng.uniformRv();

    double hastings = 0.0;
    
    if(randomMove < 0.75){
        if(!fixedTree){ // This is technically the LOCAL Move, which simultaneously updates topology and branch lengths
            moveChoice = TreeMoves::BRANCH_MOVE;
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

            std::set<Node*> neighbors1 = u->getNeighborRef();
            neighbors1.erase(v);//Exclude v
            Node* a = Node::chooseNodeFromSet(neighbors1);

            std::set<Node*> neighbors2 = v->getNeighborRef();
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
            moveChoice = TreeMoves::BRANCH_MOVE;
            branchCount += 1;
            std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
            Node* root = trees[0]->getRoot();

            Node* p = nullptr;
            do{
                p = nodes[(int)(rng.uniformRv() * nodes.size())];
            }
            while(p == root);

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
    }
    else {
        moveChoice = TreeMoves::TREE_MOVE;
        treeCount += 1;
        std::unordered_map<Node*, double> branchMapping = trees[0]->getBranchLengthMapping();
        trees[0]->updateAll();
        this->dirty();

        std::vector<double> values;
        std::vector<Node*> nodeIndices;
        double totalLength = 0.0;
        for(auto mapping : branchMapping){
            double l = mapping.second;
            nodeIndices.push_back(mapping.first);
            values.push_back(l);
            totalLength += l;
        }

        std::vector<double> alphaForward(values.size(), 0.0);
        std::vector<double> alphaReverse(values.size(), 0.0);
        std::vector<double> z(values.size(), 0.0);

        for(int i = 0; i < values.size(); i++) {
            values[i] /= totalLength;
            alphaForward[i] = (values[i] * treeAlpha) + 1;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = (z[i] * treeAlpha) + 1;
        }
        
        hastings  = Probability::Dirichlet::lnPdf(alphaReverse, values) - Probability::Dirichlet::lnPdf(alphaForward, z);

        for(int i = 0; i < values.size(); i++){
            trees[0]->setBranchLength(nodeIndices[i], z[i] * totalLength);
        }
    }

    std::vector<double> values = trees[0]->getBranchLengths();
    double totalLength = 0.0;
    for(double val : values){
        totalLength += val;
    }
    currentPrior = Probability::Gamma::lnPdf(values.size(), lambda, totalLength);

    return hastings;
}

double TreeParameter::getBranchRate() const {
    return (double)branchAcceptCount/(double)branchCount;
}

double TreeParameter::getTreeRate() const {
    return (double)treeAcceptCount/(double)treeCount;
}

void TreeParameter::tune() {
    if(branchCount > 0){
        double rate = getBranchRate();

        if ( rate > 0.33 ) {
            branchDelta *= (1.0 + ((rate-0.33)/0.67));
        }
        else {
            branchDelta /= (2.0 - rate/0.33);
        }
        branchAcceptCount = 0;
        branchCount = 0;
    }

    if(treeCount > 0){
        double rate = getTreeRate();

        if ( rate > 0.33 ) {
            treeAlpha /= (1.0 + ((rate-0.33)/0.67));
        }
        else {
            treeAlpha *= (2.0 - rate/0.33);
        }
        treeAcceptCount = 0;
        treeCount = 0;
    }
}

double TreeParameter::lnPrior() {
    return currentPrior;
}