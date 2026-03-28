#include "misc/Alignment.hpp"
#include "misc/Probability.hpp"
#include "misc/RandomVariable.hpp"
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

TreeParameter::TreeParameter(const TreeParameter& t) : priorParams{t.priorParams[0], t.priorParams[1]}, tuningState(std::make_shared<TreeTuningState>(*t.tuningState)) {

    currentPrior = t.currentPrior;
    oldPrior = t.oldPrior;
    branchAcceptCount = t.branchAcceptCount;
    branchCount = t.branchCount;
    countTuningEvents = t.countTuningEvents;
    moveChoice = t.moveChoice;
    treeAcceptCount = t.treeAcceptCount;
    treeCount = t.treeCount;

    trees[0] = new TreeObject(*t.trees[0]);
    trees[1] = new TreeObject(*t.trees[1]);

    if (t.isDirty()) {
        dirty();
    }
    else {
        clean();
    }
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
        if(countTuningEvents){
            tuningState->branchStats.acceptCount += 1;
        }
    }
    else if(moveChoice == TreeMoves::TREE_LENGTH_MOVE){
        treeAcceptCount += 1;
        if(countTuningEvents){
            tuningState->treeStats.acceptCount += 1;
        }
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
        if(countTuningEvents){
            tuningState->branchStats.count += 1;
        }
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
            alphaForward[i] = (values[i] * tuningState->branchAlpha) + 1;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = (z[i] * tuningState->branchAlpha) + 1;
        }
        
        hastings  = Probability::Dirichlet::lnPdf(alphaReverse, values) - Probability::Dirichlet::lnPdf(alphaForward, z);

        for(int i = 0; i < values.size(); i++){
            trees[0]->setBranchProportion(nodeIndices[i], z[i]);
        }
    }
    else {
        moveChoice = TreeMoves::TREE_LENGTH_MOVE;
        treeCount += 1;
        if(countTuningEvents){
            tuningState->treeStats.count += 1;
        }
        trees[0]->updateAll();

        double currentV = trees[0]->getTreeLength();
        double scale = std::exp(tuningState->treeDelta * (rng.uniformRv() - 0.5));
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
    if(tuningState->treeStats.count > 0){
        double rate = (double)tuningState->treeStats.acceptCount / (double)tuningState->treeStats.count;

        if(rate > 0.33){
            tuningState->treeDelta *= (1.0 + ((rate-0.33)/0.67));
        }
        else{
            tuningState->treeDelta /= (2.0 - rate/0.33);
        }
        tuningState->treeStats.acceptCount = 0;
        tuningState->treeStats.count = 0;
    }

    if(tuningState->branchStats.count > 0){
        double rate = (double)tuningState->branchStats.acceptCount / (double)tuningState->branchStats.count;

        if(rate > 0.33){
            tuningState->branchAlpha /= (1.0 + ((rate-0.33)/0.67));
        }
        else{
            tuningState->branchAlpha *= (2.0 - rate/0.33);
        }
        tuningState->branchStats.acceptCount = 0;
        tuningState->branchStats.count = 0;
    }
}

double TreeParameter::lnPrior() {
    return currentPrior;
}
