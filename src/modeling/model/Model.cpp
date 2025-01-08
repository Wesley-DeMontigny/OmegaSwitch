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

    int reconstructionWidth = numNodes*numChar*stateSpace;
    reconstruction = new double[reconstructionWidth];
    std::fill(reconstruction, reconstruction + reconstructionWidth, 0.0);

    int rescaleWidth = numNodes*numChar;
    rescaling = new double[rescaleWidth * 2];
    std::fill(rescaling, rescaling + rescaleWidth * 2, 0.0);

    activeT->updateAll();

    dpp->registerModel(this);
}

Model::~Model(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] reconstruction;
    delete [] activeCL;
    delete [] activeTP;
}

void Model::accept() {
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

    std::memcpy(rescaling, rescaling + numNodes*numChar, numNodes*numChar);

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
    std::vector<int> assignments = dpp->getAssignments();
    std::vector<Category> categories = dpp->getCategories();
    int numCats = dpp->getNumCategories();

    if(rateMatrix->isDirty()){
        tf::Taskflow rateTaskflow;
        activeT->updateAll();
        transProb->allocateQ(numCats);
        for(int i = 0; i < numCats; i++){
            double omega1 = categories[i].omega1;
            double omega2 = categories[i].omega2;
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
                double omega1 = categories[i].omega1;
                double omega2 = categories[i].omega2;
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
                    activeTP[nIndex] ^= true;
                    double v = activeT->getBranchLength(n);
                    for(int i = 0; i < numCats; i++)
                        transProb->setProbs(activeTP[nIndex], i, nIndex, v);
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
                double* rescalePointer = rescaling + (numChar * nIndex);
                std::fill(rescalePointer, rescalePointer + numChar, 0.0);

                for(int c = 0; c < numChar; c++){
                    double max = *pNN;
                    pNN++;
                    for(int i = 1; i < stateSpace; i++){
                        if(*pNN > max)
                            max = *pNN;
                        pNN++;
                    }
                    if(max < 1e-10){
                        pNN -= stateSpace;
                        for(int i = 1; i < stateSpace; i++){
                            *pNN /= max;
                            pNN++;
                        }
                        *rescalePointer = std::log(max);
                    }
                    rescalePointer++;
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

    double* rescalePointer = rescaling;
    for(int i = 0, len = numNodes * numChar; i < len; i++){
        lnL += *rescalePointer;
        rescalePointer++;
    }

    currentLikelihood = lnL;
}

void Model::regenerateTransitionProbs(int site, int category){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<Category> categories = dpp->getCategories();

    double omega1 = categories[category].omega1;
    double omega2 = categories[category].omega2;
    transProb->updateQ(rateMatrix->Q(omega1, omega2), category);

    for(Node* n : poSeq){
        int nIndex = n->getIndex();

        if(n != activeT->getRoot()){
            double v = activeT->getBranchLength(n);
            transProb->setProbs(activeTP[nIndex], category, nIndex, v);
        }
    }
}


double Model::testCategory(int site, int category, bool update){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<Category> categories = dpp->getCategories();

    if(update) {
        double omega1 = categories[category].omega1;
        double omega2 = categories[category].omega2;
        transProb->updateQ(rateMatrix->Q(omega1, omega2), category);
    }

    double* siteBuffer = new double[numNodes * stateSpace];
    std::fill(siteBuffer, siteBuffer + (numNodes * stateSpace), 1.0);

    double* rescaleBuffer = new double[numNodes];
    std::fill(rescaleBuffer, rescaleBuffer + numNodes, 0.0);

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

            double* rescalePointer = rescaleBuffer + nIndex;

            double max = *pNN;
            pNN++;
            for(int i = 1; i < stateSpace; i++){
                if(*pNN > max)
                    max = *pNN;
                pNN++;
            }
            if(max < 1e-10){
                pNN -= stateSpace;
                for(int i = 1; i < stateSpace; i++){
                    *pNN /= max;
                    pNN++;
                }
                *rescalePointer = std::log(max);
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
    std::vector<double> f = rateMatrix->getStationary();
    double lnL = 0.0;
    double* siteRoot = siteBuffer + (rIndex * stateSpace);

    for(int i = 0; i < stateSpace; i++){
        lnL += siteRoot[i]*f[i];
    }
    lnL = std::log(lnL);

    double* rescaleBufferPointer = rescaleBuffer;
    for(int i = 0; i < numNodes; i++){
        lnL += *rescaleBufferPointer;
        rescaleBufferPointer++;
    }

    delete [] siteBuffer;

    delete [] rescaleBuffer;

    return lnL;
}

void Model::reconstructTips(){
    TreeObject* activeT = tree->getTree();
    std::vector<Node*>&  preOrderSeq = activeT->getPostOrderSeq();
    std::reverse(preOrderSeq.begin(), preOrderSeq.end());
    std::vector<int> assignments = dpp->getAssignments();
    std::vector<Category> categories = dpp->getCategories();
    int numCats = dpp->getNumCategories();

    tf::Taskflow phyloTaskflow;
    std::unordered_map<int, tf::Task> taskMap;

    Node* root = activeT->getRoot();
    taskMap.insert(std::make_pair(root->getIndex(), phyloTaskflow.emplace([root, this, numCats, categories, assignments, activeT](){
        int rIndex = root->getIndex();
        double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
        double* rR = reconstruction + rIndex*numChar*stateSpace;
        std::vector<double> f = rateMatrix->getStationary();

        for(int c = 0; c < numChar; c++){
            double total = 0.0;
            for(int i = 0; i < stateSpace; i++){
                double product = f[i] * pR[i];
                rR[i] = product;
                total += product;
            }
            for(int i = 0; i < stateSpace; i++){
                *rR /= total;
                rR++;
            }
            pR += stateSpace;
        }
    })));

    for(Node* n : preOrderSeq){
        if(n != activeT->getRoot()){
            taskMap.insert(std::make_pair(n->getIndex(), phyloTaskflow.emplace([n, this, numCats, categories, assignments, activeT](){
                int nIndex = n->getIndex();
                double* pN = (*postOrder)(nIndex, activeCL[nIndex], 0);
                double* rN = reconstruction + nIndex*numChar*stateSpace;

                int aIndex = n->getAncestor()->getIndex();
                double* rA = reconstruction + aIndex*numChar*stateSpace;

                std::vector<Matrix<double>> transProbMatrices;
                for(int cat = 0; cat < numCats; cat++)
                    transProbMatrices.push_back(*(*transProb)(activeTP[nIndex], cat, nIndex));
                for(int c = 0; c < numChar; c++){
                    Matrix<double> P = transProbMatrices[assignments[c]];
                    double total = 0.0;
                    for(int i = 0; i < stateSpace; i++){
                        double sum = 0.0;
                        double l = pN[i];
                        for(int j = 0; j < stateSpace; j++){
                            sum += P(j, i) * rA[j] * l; //Probability of going from ancestor j to current i times probability of ancestor j times condL?
                        }
                        rN[i] = sum;
                        total += sum;
                    }

                    for(int i = 0; i < stateSpace; i++){
                        *rN /= total;
                        rN++;
                    }

                    pN += stateSpace;
                    rA += stateSpace;
                }
            })));
        }
    }

    for(Node* n : preOrderSeq){
        if(n->getIsTip() == false){
            for(Node* d : n->getNeighbors()){
                if(d != n->getAncestor()){
                    taskMap.at(n->getIndex()).precede(taskMap.at(d->getIndex()));
                }
            }
        }
    }

    executor.run(phyloTaskflow).wait();
}

void Model::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
    dpp->tune();
}

std::string Model::tabularHeader(){
    std::string returnString = "Iteration\tPosterior\tLikelihood\tTree Prior\tDPP Prior\tK Prior\tR Prior";
    returnString += "\tK\tR";
    for(int i = 0; i < 61; i++){
        returnString += "\tPi[" + std::to_string(i) + "]";
    }

    return returnString + "\n";
}

std::string Model::tabularOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                               std::to_string(currentLikelihood) + "\t" + std::to_string(tree->lnPrior()) + "\t" +
                               std::to_string(dpp->lnPrior()) + "\t" + std::to_string(rateMatrix->kPrior()) + "\t" +
                               std::to_string(rateMatrix->rPrior());
    returnString += "\t" + std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getR());
    std::vector<double> stationary = rateMatrix->getRawStationary();
    for(double i : stationary){
        returnString += "\t" + std::to_string(i);
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
        returnString += "\tOmega1[" + std::to_string(i) + "]" + "\tOmega2[" + std::to_string(i) + "]";
    return returnString + "\n";
}

std::string Model::dppOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t";
    std::vector<Category> categories = dpp->getCategories();
    returnString += std::to_string(categories.size());
    std::vector<int> assignments = dpp->getAssignments();
    for(int c : assignments){
        returnString += "\t" + std::to_string(categories[c].omega1) + "\t" + std::to_string(categories[c].omega2);
    }


    return returnString + "\n";
}

std::string Model::tipsHeader(){
    std::string returnString = "Iteration";
    std::vector<Node*> tips = tree->getTree()->getTips();

    for(Node* n : tips){
        std::string name = n->getName();
        for(int c = 0; c < numChar; c++){
            returnString += "\t" + name + "[" + std::to_string(c) + "]";
        }
    }

    return returnString + "\n";
}

std::string Model::tipsOut(int i){
    std::string returnString = std::to_string(i);
    std::vector<Node*> tips = tree->getTree()->getTips();

    std::vector<int> assignments = dpp->getAssignments();
    std::vector<Category> categories = dpp->getCategories();

    for(Node* n : tips) {
        double* rN = reconstruction + n->getIndex()*numChar*stateSpace;
        for(int c = 0; c < numChar; c++) {
            double expectedOmega = 0;
            double omega1 = categories[assignments[c]].omega1;
            double omega2 = categories[assignments[c]].omega2;

            for(int i = 0; i < stateSpace; i++) {
                if(*rN > 0) {
                    if(stateSpace < 61) {
                        expectedOmega += omega1 * (*rN);
                    }
                    else {
                        expectedOmega += omega2 * (*rN);
                    }
                }

                rN++;
            }

            returnString += "\t" + std::to_string(expectedOmega);
        }
    }

    return returnString + "\n";
}