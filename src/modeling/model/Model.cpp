#include "Model.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/CodonMultiMatrix.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>

Model::Model(Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d) : aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), dpp(d) {

    TreeObject* activeT = tree->getTree();
    stateSpace = 122;
    int numChar = aln->getNumChar();

    dpp->registerModel(this);

    if(aln->getNumTaxa() != activeT->getNumTaxa())
        Msg::error("Expected " + std::to_string(aln->getNumTaxa()) + 
        "taxa in the tree, but found only " + std::to_string(activeT->getNumTaxa()));
    
    int numNodes = (aln->getNumTaxa() * 2) - 1;

    int flagWidths = 2 * numNodes;
    activeCL = new bool[flagWidths];
    activeTP = new bool[flagWidths];
    for(int i = 0; i < flagWidths; i++){
        activeCL[i] = false;
        activeTP[i] = false;
    }

    //We need to do some setting to make sure the alignment and tree match
    std::vector<std::string> taxaNames = aln->getTaxaNames();
    bool randomAssign = false;
    for(Node* n : activeT->getTips()){
        bool found = false;
        std::string name = n->getName();
        for(int i = 0; i < taxaNames.size(); i++){
            if(name == taxaNames[i]){
                found=true;
                n->setIndex(i);
                break;
            }
        }

        if(found == false){
            Msg::warning("Expected to find a sequence named " + name + "! Assigning indices and names randomly.");
            randomAssign = true;
            break;
        }
    }

    //Set the tips randomly if the tips are not properly named
    if(randomAssign == true){
        std::vector<Node*> tips = activeT->getTips();
        for(int i = 0; i < taxaNames.size(); i++){
            tips[i]->setIndex(i);
            tips[i]->setName(taxaNames[i]);
        }
    }

    tree->accept(); //Accept the tip changes into memory tree (if any happened)

    postOrder = new ConditionalLikelihood(aln, 1);
    transProb = new TransitionProbability(numNodes, numChar);

    activeT->updateAll();
}

Model::~Model(){
    delete postOrder;
    delete transProb;
    delete [] activeCL;
    delete [] activeTP;
}

void Model::accept() {
    oldLikelihood = currentLikelihood;

    int numNodes = (aln->getNumTaxa() * 2) - 1;
    for(int i = 0; i < numNodes; i++){
        activeCL[i + numNodes] = activeCL[i];
        activeTP[i + numNodes] = activeTP[i];
    }

    if(tree->isDirty()){
        tree->accept();
        tree->clean();
    }
    if(rateMatrix->isDirty()){
        rateMatrix->accept();
        rateMatrix->clean();
    }
    if(dpp->isDirty()){
        dpp->accept();
        dpp->clean();
    }

    transProb->accept();
}

void Model::reject() {
    currentLikelihood = oldLikelihood;

    int numNodes = (aln->getNumTaxa() * 2) - 1;
    for(int i = 0; i < numNodes; i++){
        activeCL[i] = activeCL[i + numNodes];
        activeTP[i] = activeTP[i + numNodes];
    }

    if(tree->isDirty()){
        tree->reject();
        tree->clean();
    }
    if(rateMatrix->isDirty()){
        rateMatrix->reject();
        rateMatrix->clean();
    }
    if(dpp->isDirty()){
        dpp->reject();
        dpp->clean();
    }

    transProb->reject();
}

double Model::lnPrior(){
    return dpp->lnPrior() + tree->lnPrior() + rateMatrix->lnPrior();
}

void Model::postOrderPrune(){

    TreeObject* activeT = tree->getTree();
    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<int> assignments = dpp->getAssinments();

    int numCats = dpp->getNumCategories();
    int numChar = aln->getNumChar();

    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        //Only update the conditional likelihoods if the node has changed
        if(n->getNeedsTPUpdate() == true){
            if(n != activeT->getRoot()){
                activeTP[nIndex] ^= true;
                double v = activeT->getBranchLength(n);
                for(int i = 0; i < numCats; i++){
                    transProb->setProbs(activeTP[nIndex], i, nIndex, v);
                }
            }
            n->setNeedsTPUpdate(false);
        }
        if(n->getNeedsCLUpdate() == true || n->getNeedsTPUpdate() == true){
            //Get memory address of the node we are looking at and pre-set all of the likelihoods at each site to be 1.0
            activeCL[nIndex] ^= true;
            double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0);
            std::fill(pNN, pNN + (numChar* stateSpace), 1.0);

            std::set<Node*>& nNeighbors = n->getNeighbors();
            //Iterate over the descendents (usually only two)
            for(Node* d : nNeighbors){
                if(d != n->getAncestor()){
                    int dIndex = d->getIndex();
                    double* pN = pNN;
                    double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0);

                    //Iterate over each of the characters and each of the potential states of our node
                    for(int c = 0; c < numChar; c++){
                        Matrix<double> P = *(*transProb)(activeTP[dIndex], assignments[c], dIndex);
                        for(int i = 0; i < stateSpace; i++){
                            //Sum up the products of the likelihoods from the CTMC (transitioning from the node's hypothetical state to another) and the conditional likelihood of the descendent states 
                            double sum = 0.0;
                            for(int j = 0; j < stateSpace; j++){
                                sum += P(i, j) * pD[j];
                            }
                            //If this is the first time, set pN equal to the sum, otherwise multiply them
                            (*pN) *= sum;
                            //Move the memory address to the next character state
                            pN++;
                        }
                        //Move the memory address to the next site
                        pD+=stateSpace;
                    }
                }
            }
            //Note that we have updated the node
            n->setNeedsCLUpdate(false);
        }
    }
}

void Model::regenerateLikelihood(){
    TreeObject* activeT = tree->getTree();
    int numCats = dpp->getNumCategories();

    if(rateMatrix->isDirty() || dpp->isDirty()){
        activeT->updateAll();
        for(int i = 0; i < numCats; i++){
            double omega1 = dpp->getCategoryOmega(i);
            double omega2 = omega1 * dpp->getCategoryBeta(i);
            transProb->updateQ(rateMatrix->Q(omega1, omega2), i);
        }
    }
    postOrderPrune();

    //Calculate the likelihood of the tree by summing up the likelihood at the root.
    int numChar = aln->getNumChar();
    int rIndex = activeT->getRoot()->getIndex();
    double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
    std::vector<double> f = rateMatrix->getStationary();
    double lnL = 0.0;

    for(int c = 0; c < numChar; c++){
        double like = 0.0;
        for(int i = 0; i < stateSpace; i++){
            like += pR[i]*f[i];
        }

        lnL += std::log(like);
        pR += stateSpace;
    }

    currentLikelihood = lnL;
}

// ONLY DO THIS IF YOU CAN PROMISE IT WILL BE ACCEPTED!!
void Model::forceRegenerate(int site, int category, bool update) {
    int numChar = aln->getNumChar();

    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    if(update){
        double omega1 = dpp->getCategoryOmega(category);
        double omega2 = omega1 * dpp->getCategoryBeta(category);
        transProb->updateQ(rateMatrix->Q(omega1, omega2), category);
    }

    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        
        if(n != activeT->getRoot()){
            double v = activeT->getBranchLength(n);
            transProb->setProbs(activeTP[nIndex], category, nIndex, v);
        }

        double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0) + (site*stateSpace);
        std::fill(pNN, pNN + stateSpace, 1.0);

        std::set<Node*>& nNeighbors = n->getNeighbors();

        for(Node* d : nNeighbors){
            if(d != n->getAncestor()){
                int dIndex = d->getIndex();
                double* pN = pNN;
                double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0) + (site*stateSpace);
            
                //Regenerate site specific likelihood
                Matrix<double> P = *(*transProb)(activeTP[dIndex], category, dIndex);
                for(int i = 0; i < stateSpace; i++){
                    //Sum up the products of the likelihoods from the CTMC (transitioning from the node's hypothetical state to another) and the conditional likelihood of the descendent states 
                    double sum = 0.0;
                    for(int j = 0; j < stateSpace; j++){
                        sum += P(i, j) * pD[j];
                    }
                    //If this is the first time, set pN equal to the sum, otherwise multiply them
                    (*pN) *= sum;
                    //Move the memory address to the next character state
                    pN++;
                }
            }
        }
    }

    //Calculate the likelihood of the tree by summing up the likelihood at the root.
    int rIndex = activeT->getRoot()->getIndex();
    double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
    double* pRSite = (*postOrder)(rIndex, activeCL[rIndex], 0) + (site * stateSpace);
    std::vector<double> f = rateMatrix->getStationary();
    double lnL = 0.0;

    for(int c = 0; c < numChar; c++){
        if(c != site){
            double like = 0.0;
            for(int i = 0; i < stateSpace; i++){
                like += pR[i]*f[i];
            }
            lnL += std::log(like);
            pR += stateSpace;
        }
        else{
            double like = 0.0;
            for(int i = 0; i < stateSpace; i++){
                like += pRSite[i]*f[i];
            }
            lnL += std::log(like);
        }
    }
}

double Model::regenerateIntoLikelihoodBuffer(int site, int category, bool update){

    int numChar = aln->getNumChar();

    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    if(update){
        double omega1 = dpp->getCategoryOmega(category);
        double omega2 = omega1 * dpp->getCategoryBeta(category);
        transProb->updateQ(rateMatrix->Q(omega1, omega2), category);
    }

    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        
        if(n != activeT->getRoot()){
            double v = activeT->getBranchLength(n);
            transProb->setProbs(activeTP[nIndex], category, nIndex, v);
        }

        double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0) + (numChar*stateSpace);
        std::fill(pNN, pNN + stateSpace, 1.0);

        std::set<Node*>& nNeighbors = n->getNeighbors();

        for(Node* d : nNeighbors){
            if(d != n->getAncestor()){
                int dIndex = d->getIndex();
                double* pN = pNN;
                double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0) + (numChar*stateSpace);
            
                //Regenerate site specific likelihood
                Matrix<double> P = *(*transProb)(activeTP[dIndex], category, dIndex);
                for(int i = 0; i < stateSpace; i++){
                    //Sum up the products of the likelihoods from the CTMC (transitioning from the node's hypothetical state to another) and the conditional likelihood of the descendent states 
                    double sum = 0.0;
                    for(int j = 0; j < stateSpace; j++){
                        sum += P(i, j) * pD[j];
                    }
                    //If this is the first time, set pN equal to the sum, otherwise multiply them
                    (*pN) *= sum;
                    //Move the memory address to the next character state
                    pN++;
                }
            }
        }
    }

    //Calculate the likelihood of the tree by summing up the likelihood at the root.
    int rIndex = activeT->getRoot()->getIndex();
    double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
    double* pRSite = (*postOrder)(rIndex, activeCL[rIndex], 0) + (numChar * stateSpace);
    std::vector<double> f = rateMatrix->getStationary();
    double lnL = 0.0;

    for(int c = 0; c < numChar; c++){
        if(c != site){
            double like = 0.0;
            for(int i = 0; i < stateSpace; i++){
                like += pR[i]*f[i];
            }
            lnL += std::log(like);
            pR += stateSpace;
        }
        else{
            double like = 0.0;
            for(int i = 0; i < stateSpace; i++){
                like += pRSite[i]*f[i];
            }
            lnL += std::log(like);
        }
    }

    return lnL;
}

std::string Model::tabularOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                               std::to_string(currentLikelihood) + "\t" + std::to_string(tree->lnPrior()) + "\t" +
                               std::to_string(dpp->lnPrior()) + "\t" + std::to_string(rateMatrix->kPrior()) + "\t" +
                               std::to_string(rateMatrix->rPrior());
    if(rateMatrix->updatingStationary())
        returnString += "\t" + std::to_string(rateMatrix->stationaryPrior());
    returnString += "\t" + std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getR());
    if(rateMatrix->updatingStationary()){
        std::vector<double> stationary = rateMatrix->getStationary();
        for(double i : stationary){
            returnString += "\t" + std::to_string(i);
        }
    }

    return returnString;
}

std::string Model::tabularHeader(){
    std::string returnString = "Iteration\tPosterior\tLikelihood\tTree Prior\tDPP Prior\tK Prior\tR Prior";
    if(rateMatrix->updatingStationary())
        returnString += "\tPi Prior";
    returnString += "\tK\tR";
    if(rateMatrix->updatingStationary()){
        for(int i = 0; i < 122; i++){
            returnString += "\tPi[" + std::to_string(i) + "]";
        }
    }

    return returnString;
}

std::string Model::treeOut(int i){
    return "";
}

std::string Model::treeHeader(){
    return "";
}

std::string Model::dppOut(int i){
    return "";
}

std::string Model::dppHeader(){
    return "";
}