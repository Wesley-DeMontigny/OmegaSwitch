#include "TreeParameter.hpp"
#include "core/RandomVariable.hpp"
#include "TreeObject.hpp"

TreeParameter::TreeParameter(int nt){
    trees[0] = new TreeObject(nt);
    trees[1] = new TreeObject(*trees[0]);
    dirty();
}

TreeParameter::TreeParameter(Alignment* aln){
    trees[0] = new TreeObject(aln);
    trees[1] = new TreeObject(*trees[0]);
    dirty();
}

TreeParameter::TreeParameter(std::string newick, std::vector<std::string> taxaNames){
    trees[0] = new TreeObject(newick, taxaNames);
    trees[1] = new TreeObject(*trees[0]);
    dirty();
}

TreeParameter::~TreeParameter(){
    delete trees[0];
    delete trees[1]; 
}

TreeParameter& TreeParameter::operator=(const TreeParameter& t){
    if(this == &t)
        return *this;

    const TreeObject* activeTree = t.getTreeConst();
    *trees[0] = *activeTree;
    *trees[1] = *activeTree;

    return *this;
}

void TreeParameter::accept(){
    *trees[1] = *trees[0];
}

std::vector<double> TreeParameter::getBranchLengths(){
    return trees[0]->getBranchLengths();
}

void TreeParameter::reject(){
    *trees[0] = *trees[1];
}

void TreeParameter::regenerate() {}

std::string TreeParameter::writeValue(){
    return trees[0]->getNewick();
}