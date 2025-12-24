#include "core/Alignment.hpp"
#include "core/Probability.hpp"
#include "core/RandomVariable.hpp"
#include "Node.hpp"
#include "TreeObject.hpp"
#include "TreeParameter.hpp"
#include <cmath>

TreeParameter::TreeParameter(Alignment& aln, std::string& newick, double p[2]) : priorParams{p[0], p[1]} {

    trees[0] = new TreeObject(newick, aln.getTaxaNames());

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> nodes = trees[0]->getPostOrderSeq();

    double treeLength = Probability::Gamma::rv(&rng, priorParams[0], priorParams[1]);
    trees[0]->setTreeLength(treeLength);
    currentPrior = Probability::Gamma::lnPdf(priorParams[0], priorParams[1], treeLength);
    oldPrior = currentPrior;

    trees[1] = new TreeObject(*trees[0]);

    dirty();
}

TreeParameter::TreeParameter(TreeObject& tree, double p[2]) : priorParams{p[0], p[1]} {

    trees[0] = new TreeObject(tree);

    // Randomize these branch lengths
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    trees[0]->randomizeBranches();
    double treeLength = Probability::Gamma::rv(&rng, priorParams[0], priorParams[1]);
    trees[0]->setTreeLength(treeLength);
    currentPrior = Probability::Gamma::lnPdf(priorParams[0], priorParams[1], treeLength);
    oldPrior = currentPrior;

    trees[1] = new TreeObject(*trees[0]);

    dirty();
}

TreeParameter::~TreeParameter(){
    delete trees[0];
    delete trees[1]; 
}

void TreeParameter::accept(){
    *trees[1] = *trees[0];
    oldPrior = currentPrior;

    if(moveChoice == TreeMoves::BRANCH_PROPORTION_MOVE){
        branchAcceptCount += 1;
    }
    else if(moveChoice == TreeMoves::TREE_LENGTH_MOVE){
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
    this->dirty();
    
    if(randomMove < 0.75){
        moveChoice = TreeMoves::BRANCH_PROPORTION_MOVE;
        branchCount += 1;
        std::unordered_map<Node*, double> branchMapping = trees[0]->getBranchPropMapping();
        trees[0]->updateAll();

        std::vector<double> values;
        std::vector<Node*> nodeIndices;
        for(auto mapping : branchMapping){
            double l = mapping.second;
            nodeIndices.push_back(mapping.first);
            values.push_back(l);
        }

        std::vector<double> alphaForward(values.size(), 0.0);
        std::vector<double> alphaReverse(values.size(), 0.0);
        std::vector<double> z(values.size(), 0.0);

        for(int i = 0; i < values.size(); i++) {
            alphaForward[i] = (values[i] * branchAlpha) + 1;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = (z[i] * branchAlpha) + 1;
        }
        
        hastings  = Probability::Dirichlet::lnPdf(alphaReverse, values) - Probability::Dirichlet::lnPdf(alphaForward, z);

        for(int i = 0; i < values.size(); i++){
            trees[0]->setBranchProportion(nodeIndices[i], z[i]);
        }
    }
    else {
        moveChoice = TreeMoves::TREE_LENGTH_MOVE;
        treeCount += 1;
        trees[0]->updateAll();

        double currentV = trees[0]->getTreeLength();
        double scale = std::exp(treeDelta * (rng.uniformRv() - 0.5));
        double newV = currentV * scale;

        trees[0]->setTreeLength(newV);
        hastings = std::log(scale);

        currentPrior = Probability::Gamma::lnPdf(priorParams[0], priorParams[1], newV);
    }

    return hastings;
}

double TreeParameter::getBranchRate() const {
    return (double)branchAcceptCount/(double)branchCount;
}

double TreeParameter::getTreeRate() const {
    return (double)treeAcceptCount/(double)treeCount;
}

void TreeParameter::tune() {
    if(treeCount > 0){
        double rate = getTreeRate();

        if ( rate > 0.33 ) {
            treeDelta *= (1.0 + ((rate-0.33)/0.67));
        }
        else {
            treeDelta /= (2.0 - rate/0.33);
        }
        treeAcceptCount = 0;
        treeCount = 0;
    }

    if(branchCount > 0){
        double rate = getBranchRate();

        if ( rate > 0.33 ) {
            branchAlpha /= (1.0 + ((rate-0.33)/0.67));
        }
        else {
            branchAlpha *= (2.0 - rate/0.33);
        }
        branchAcceptCount = 0;
        branchCount = 0;
    }
}

double TreeParameter::lnPrior() {
    return currentPrior;
}