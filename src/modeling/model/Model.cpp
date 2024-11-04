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
#include <unordered_map>
#include <chrono>

Model::Model(Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d) : aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), dpp(d), numChar(0) {

    TreeObject* activeT = tree->getTree();
    stateSpace = 122;
    numChar = aln->getNumChar();

    dpp->registerModel(this);

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
    postOrder = new ConditionalLikelihood(aln, numNodes, 1);
    transProb = new TransitionProbability(numNodes, numChar + 5);

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

void Model::regenerateLikelihood(){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<int> assignments = dpp->getAssinments();
    std::vector<Category> categories = dpp->getCategories();
    int numCats = dpp->getNumCategories();

    if(rateMatrix->isDirty()){
        tf::Taskflow rateTaskflow;
        activeT->updateAll();
        transProb->allocateQ(numCats);
        for(int i = 0; i < numCats; i++){
            double omega1 = categories[i].omega;
            double omega2 = omega1 * categories[i].beta;
            rateTaskflow.emplace([this, omega1, omega2, i](){
                transProb->updateQ(rateMatrix->Q(omega1, omega2), i);
            });
        }

        executor.run(rateTaskflow).wait();
    }
    else if(dpp->isDirty()){
        activeT->updateAll();
        transProb->allocateQ(numCats);
        for(int i = 0; i < numCats; i++){
            if(categories[i].dirty){
                double omega1 = categories[i].omega;
                double omega2 = omega1 * categories[i].beta;
                transProb->updateQ(rateMatrix->Q(omega1, omega2), i);
            }
        }
    }

    tf::Taskflow phyloTaskflow;
    std::unordered_map<int, tf::Task> taskMap;

    for(Node* n : poSeq){
        taskMap.insert(std::make_pair(n->getIndex(), phyloTaskflow.emplace([n, this, numCats, categories, assignments, activeT](){
            int nIndex = n->getIndex();
            if(n->getNeedsTPUpdate() == true){
                if(n != activeT->getRoot()) {
                    if(!dpp->isDirty()) {
                        activeTP[nIndex] ^= true;
                        double v = activeT->getBranchLength(n);
                        for(int i = 0; i < numCats; i++)
                            transProb->setProbs(activeTP[nIndex], i, nIndex, v);
                    }
                    else {
                        activeTP[nIndex] ^= true;
                        double v = activeT->getBranchLength(n);
                        for(int i = 0; i < numCats; i++){
                            if(categories[i].dirty)
                                transProb->setProbs(activeTP[nIndex], i, nIndex, v);
                            else
                                transProb->pullProbs(activeTP[nIndex], i, nIndex, v);
                        }
                    }
                }
                n->setNeedsTPUpdate(false);
            }
            if(n->getNeedsCLUpdate() == true){
                activeCL[nIndex] ^= true;
                double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0);
                std::fill(pNN, pNN + (numChar * stateSpace), 1.0);

                std::set<Node*>& nNeighbors = n->getNeighbors();
                for(Node* d : nNeighbors){
                    if(d != n->getAncestor()){
                        int dIndex = d->getIndex();
                        double* pN = pNN;
                        double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0);

                        std::vector<Matrix<double>> transProbMatrices;
                        for(int cat = 0; cat < numCats; cat++)
                            transProbMatrices.push_back(*(*transProb)(activeTP[dIndex], cat, dIndex));
                        for(int c = 0; c < numChar; c++){
                            Matrix<double> P = transProbMatrices[assignments[c]];
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
                n->setNeedsCLUpdate(false);
            }
        })));
    }

    for(Node* n : poSeq){
        if(n != activeT->getRoot())
            taskMap.at(n->getIndex()).precede(taskMap.at(n->getAncestor()->getIndex()));
    }

    executor.run(phyloTaskflow).wait();

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

void Model::regenerateLikelihood(int site, int category, bool update){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<Category> categories = dpp->getCategories();

    if(update){
        double omega1 = categories[category].omega;
        double omega2 = omega1 * categories[category].beta;
        transProb->updateQ(rateMatrix->Q(omega1, omega2), category);
    }

    for(Node* n : poSeq){
        int nIndex = n->getIndex();

        if(update){
            if(n != activeT->getRoot()){
                double v = activeT->getBranchLength(n);
                transProb->setProbs(activeTP[nIndex], category, nIndex, v);
            }
        }

        if(n->getIsTip() == false){
            double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0) + site*stateSpace;
            std::fill(pNN, pNN + stateSpace, 1.0);

            std::set<Node*>& nNeighbors = n->getNeighbors();

            for(Node* d : nNeighbors){
                if(d != n->getAncestor()){
                    int dIndex = d->getIndex();
                    double* pN = pNN;
                    double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0) + site*stateSpace;

                    Matrix<double> P = *(*transProb)(activeTP[dIndex], category, dIndex);
                    for(int i = 0; i < stateSpace; i++){
                        double sum = 0.0;
                        for(int j = 0; j < stateSpace; j++){
                            sum += P(i, j) * pD[j];
                        }
                        (*pN) *= sum;
                        pN++;
                    }
                }
            }
        }
    }

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


double Model::testCategory(int site, int category, bool update){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<Category> categories = dpp->getCategories();

    if(update){
        double omega1 = categories[category].omega;
        double omega2 = omega1 * categories[category].beta;
        transProb->updateQ(rateMatrix->Q(omega1, omega2), category);
    }

    double* siteBuffer = new double[numNodes * stateSpace];
    std::fill(siteBuffer, siteBuffer + (numNodes * stateSpace), 1.0);

    for(Node* n : poSeq){
        int nIndex = n->getIndex();

        if(update){
            if(n != activeT->getRoot()){
                double v = activeT->getBranchLength(n);
                transProb->setProbs(activeTP[nIndex], category, nIndex, v);
            }
        }
        
        double* pNN = siteBuffer + (nIndex * stateSpace);
        if(n->getIsTip() == false){

            std::set<Node*>& nNeighbors = n->getNeighbors();

            for(Node* d : nNeighbors){
                if(d != n->getAncestor()){
                    int dIndex = d->getIndex();
                    double* pN = pNN;
                    double* pD = siteBuffer + (dIndex * stateSpace);

                    Matrix<double> P = *(*transProb)(activeTP[dIndex], category, dIndex);
                    for(int i = 0; i < stateSpace; i++){
                        double sum = 0.0;
                        for(int j = 0; j < stateSpace; j++){
                            sum += P(i, j) * pD[j];
                        }
                        (*pN) *= sum;
                        pN++;
                    }
                }
            }
        }
        else {
            double* pTip = (*postOrder)(nIndex, activeCL[nIndex], 0) + site*stateSpace;
            for(int i = 0; i < stateSpace; i++){
                *pNN = *pTip;
                pTip++;
                pNN++;
            }
        }
    }

    int rIndex = activeT->getRoot()->getIndex();
    double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
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
        else {
            double like = 0.0;
            double* siteRoot = siteBuffer + (rIndex * stateSpace);
            for(int i = 0; i < stateSpace; i++){
                like += siteRoot[i]*f[i];
            }
            lnL += std::log(like);
            pR += stateSpace;
        }
    }

    delete [] siteBuffer;

    return lnL;
}


void Model::tuneMoves(){
    dpp->tune();
    tree->tune();
    rateMatrix->tune();
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

    return returnString + "\n";
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

    return returnString + "\n";
}

std::string Model::treeHeader(){
    return "Iteration\tPosterior\tTree\n";
}

std::string Model::treeOut(int i){
    return std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" + tree->writeNewick() + "\n";
}

std::string Model::dppHeader(){
    std::string returnString = "Iteration\tPosterior\tCategoryCount";
    for(int i = 0; i < numChar; i++)
        returnString += "\tOmega[" + std::to_string(i) + "]" + "\tBeta[" + std::to_string(i) + "]";

    return returnString + "\n";
}

std::string Model::dppOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t";
    std::vector<Category> categories = dpp->getCategories();
    returnString += std::to_string(categories.size());
    std::vector<int> assignments = dpp->getAssinments();
    for(int c : assignments){
        returnString += "\t" + std::to_string(categories[c].omega) + "\t" + std::to_string(categories[c].beta);
    }

    return returnString + "\n";
}
