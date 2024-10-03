#include "PhyloCTMC.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/rates/CodonMultiMatrix.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>

PhyloCTMC::PhyloCTMC(Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d) : aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), dpp(d) {

    this->dirty();

    TreeObject* activeT = tree->getTree();
    stateSpace = 122;
    if(aln->getNumTaxa() != activeT->getNumTaxa())
        Msg::error("Expected " + std::to_string(aln->getNumTaxa()) + 
        "taxa in the tree, but found only " + std::to_string(activeT->getNumTaxa()));
    
    int numNodes = aln->getNumTaxa() * 2;
    activeCL = new bool[2 * numNodes];
    activeTP = new bool[2 * numNodes];
    for(int i = 0; i < 2 * numNodes; i++){
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
    transProb = new TransitionProbability(numNodes, aln->getNumChar());

    activeT->updateAll();
}

PhyloCTMC::~PhyloCTMC(){
    delete postOrder;
    delete transProb;
    delete [] activeCL;
    delete [] activeTP;
}

void PhyloCTMC::accept() {
    oldLikelihood = currentLikelihood;

    int numNodes = aln->getNumTaxa() * 2;
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

    transProb->accept();
}

void PhyloCTMC::reject() {
    currentLikelihood = oldLikelihood;

    int numNodes = aln->getNumTaxa() * 2;
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

    transProb->reject();
}

void PhyloCTMC::regenerate(){
    tree->regenerate();
    rateMatrix->regenerate();

    if(tree->isDirty() || rateMatrix->isDirty())
        this->dirty();

    int numChar = aln->getNumChar();

    if(this->isDirty()){
        TreeObject* activeT = tree->getTree();

        std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
        std::vector<int> assignments = dpp->getAssinments();
        int numCats = dpp->getNumCategories();

        if(rateMatrix->isDirty()){
            activeT->updateAll();
            for(int i = 0; i < numCats; i++){
                transProb->updateQ(rateMatrix->Q(i), i);
            }
        }

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

        //Calculate the likelihood of the tree by summing up the likelihood at the root.
        int rIndex = activeT->getRoot()->getIndex();
        double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
        std::vector<double> f = rateMatrix->stationary();
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
}

double PhyloCTMC::regenerateAtSite(int site, int category, bool update){
    this->dirty();

    int numChar = aln->getNumChar();

    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    if(update)
        transProb->updateQ(rateMatrix->Q(category), category);

    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        
        if(n != activeT->getRoot()){
            double v = activeT->getBranchLength(n);
            transProb->setProbs(activeTP[nIndex] ^ true, category, nIndex, v);
        }

        double* pNN = (*postOrder)(nIndex, activeCL[nIndex] ^ true, 0) + (site*stateSpace);
        std::fill(pNN, pNN + stateSpace, 1.0);

        std::set<Node*>& nNeighbors = n->getNeighbors();

        for(Node* d : nNeighbors){
            if(d != n->getAncestor()){
                int dIndex = d->getIndex();
                double* pN = pNN;
                double* pD = (*postOrder)(dIndex, activeCL[dIndex] ^ true, 0) + (site*stateSpace);
            
                //Regenerate site specific likelihood
                Matrix<double> P = *(*transProb)(activeTP[dIndex] ^ true, category, dIndex);
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
    double* pRSite = (*postOrder)(rIndex, activeCL[rIndex] ^ true, 0) + (site * stateSpace);
    std::vector<double> f = rateMatrix->stationary();
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