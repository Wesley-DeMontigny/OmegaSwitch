#include "M0Model.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/M0Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <unordered_map>
//#include <chrono>

M0Model::M0Model(Settings s, Alignment* a, TreeParameter* t, M0Matrix* m) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0),numChar(0) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    TreeObject* activeT = tree->getTree();
    stateSpace = 61;
    numChar = aln->getNumChar();

    if(aln->getNumTaxa() != activeT->getNumTaxa())
        Msg::error("Expected " + std::to_string(aln->getNumTaxa()) + 
        "taxa in the tree, but found only " + std::to_string(activeT->getNumTaxa()));
    
    numNodes = tree->getTree()->getNumNodes();

    int flagWidths = 2 * numNodes;
    activeCL = new bool[flagWidths];
    activeTP = new bool[flagWidths];
    for(int i = 0; i < flagWidths; i++){
        activeCL[i] = false;
        activeTP[i] = false;
    }
    postOrder = new ConditionalLikelihood(aln, numNodes, 1, stateSpace);
    transProb = new TransitionProbability(numNodes, stateSpace);

    int rescaleWidth = numNodes*numChar;
    rescaling = new double[rescaleWidth * 2];
    std::fill(rescaling, rescaling + rescaleWidth * 2, 0.0);

    activeT->updateAll();
}

M0Model::~M0Model(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] activeCL;
    delete [] activeTP;
}

void M0Model::accept() {
    oldLikelihood = currentLikelihood;

    for(int i = 0; i < numNodes; i++){
        activeCL[i + numNodes] = activeCL[i];
        activeTP[i + numNodes] = activeTP[i];
    }

    std::memcpy(rescaling + numNodes*numChar, rescaling, numNodes*numChar);

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

void M0Model::reject() {
    currentLikelihood = oldLikelihood;

    for(int i = 0; i < numNodes; i++){
        activeCL[i] = activeCL[i + numNodes];
        activeTP[i] = activeTP[i + numNodes];
    }

    std::memcpy(rescaling, rescaling + numNodes*numChar, numNodes*numChar);

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

double M0Model::lnPrior(){
    return tree->lnPrior() + rateMatrix->lnPrior();
}

void M0Model::regenerateLikelihood(){
    TreeObject* activeT = tree->getTree();

    const std::vector<Node*> poSeq = activeT->getPostOrderSeq();

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if(rateMatrix->isDirty()){
        activeT->updateAll();
        transProb->updateQ(rateMatrix->Q(), 0);
    }

    //std::chrono::steady_clock::time_point rateTime = std::chrono::steady_clock::now();
    //std::cout << "Rate computation was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(rateTime - begin).count() << "[milliseconds]" << std::endl;

    tf::Taskflow probsTaskflow;

    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        if(n->getNeedsTPUpdate() == true){
            if(n != activeT->getRoot()) {
                double v = activeT->getBranchLength(n);
                activeTP[nIndex] ^= true;
                bool activeIndex = activeTP[nIndex];
                probsTaskflow.emplace([this, nIndex, v, activeIndex](){
                    transProb->setProbs(activeIndex, 0, nIndex, v);
                });
            }
            n->setNeedsTPUpdate(false);
        }
    }

    executor.run(probsTaskflow).wait();

    //std::chrono::steady_clock::time_point probsTime = std::chrono::steady_clock::now();
    //std::cout << "Probs computation was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(probsTime - rateTime).count() << "[milliseconds]" << std::endl;

    for(Node* n : poSeq){
        if(n->getNeedsCLUpdate() == true){
            activeCL[n->getIndex()] ^= true; // Flip this ahead of time
        }
    }

    tf::Taskflow phyloTaskflow;
    
    int chunkSize = 50;
    for(int range = 0; range < (int)std::ceil((double)numChar / chunkSize); range++){
        int start = range * chunkSize;
        int end = start + chunkSize-1;
        end = std::min(end, numChar-1);

        phyloTaskflow.emplace([this, &poSeq, start, end](){
            int currentChunkSize = end - start + 1;
            for(Node* n : poSeq){
                int nIndex = n->getIndex();
                if(n->getNeedsCLUpdate() == true){
                    double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0) + start * stateSpace;
                    std::fill(pNN, pNN + (currentChunkSize * stateSpace), 1.0);

                    std::set<Node*>& nNeighbors = n->getNeighbors();
                    for(Node* d : nNeighbors){
                        if(d != n->getAncestor()){
                            int dIndex = d->getIndex();
                            double* pN = pNN;
                            double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0) + start * stateSpace;

                            const Matrix<double>& P = (*transProb)(activeTP[dIndex], 0, dIndex);
                            for(int c = 0; c < currentChunkSize; c++){
                                for(int i = 0; i < stateSpace; i++){
                                    double sum = 0.0;
                                    for(int j = 0; j < stateSpace; j++){
                                        sum += P(i, j) * pD[j];
                                    }
                                    (*pN) *= sum;
                                    pN++;
                                }
                                pD+=stateSpace;
                            }
                        }
                    }
                    double* rescalePointer = rescaling + (numChar * nIndex) + start;
                    std::fill(rescalePointer, rescalePointer + currentChunkSize, 0.0);

                    for(int c = 0; c < currentChunkSize; c++){
                        double max = *pNN;
                        pNN++;
                        for(int i = 1; i < stateSpace; i++){
                            if(*pNN > max)
                                max = *pNN;
                            pNN++;
                        }
                        if(max < 1e-10){
                            pNN -= stateSpace;
                            for(int i = 0; i < stateSpace; i++){
                                *pNN /= max;
                                pNN++;
                            }
                            *rescalePointer = std::log(max);
                        }
                        rescalePointer++;
                    }
                }
            }
        });
    }
    executor.run(phyloTaskflow).wait();

    // Mark everything updated
    for(Node* n : poSeq)
        n->setNeedsCLUpdate(false);

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

    double* rescalePointer = rescaling;
    for(int i = 0, len = numNodes * numChar; i < len; i++){
        lnL += *rescalePointer;
        rescalePointer++;
    }

    currentLikelihood = lnL;

    //std::chrono::steady_clock::time_point pruneTime = std::chrono::steady_clock::now();
    //std::cout << "Pruning was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(pruneTime - probsTime).count() << "[milliseconds]" << std::endl;
}

void M0Model::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
}

std::string M0Model::tabularHeader(){
    std::string returnString = "Iteration\tPosterior\tLikelihood\tPrior\tK\tR\tOmega";
    for(int i = 0; i < 61; i++){
        returnString += "\tPi[" + std::to_string(i) + "]";
    }

    return returnString + "\n";
}

std::string M0Model::tabularOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                               std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" +
                               std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getOmega());
    std::vector<double> stationary = rateMatrix->getRawStationary();
    for(double i : stationary){
        returnString += "\t" + std::to_string(i);
    }

    return returnString + "\n";
}

std::string M0Model::treeHeader(){
    return "Iteration\tPosterior\tTree\n";
}

std::string M0Model::treeOut(int i){
    return std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" + tree->writeNewick() + "\n";
}