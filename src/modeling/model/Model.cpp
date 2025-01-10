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
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <unordered_map>
#include <chrono>

Model::Model(Settings s, Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0),
            dpp(d), numChar(0), invariantUpdate(false), numGibbsUpdate(s.numGibbsUpdate) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

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
    transProb = new TransitionProbability(numNodes, (numChar + 5)*2);

    int rescaleWidth = numNodes*numChar;
    rescaling = new double[rescaleWidth * 2];
    std::fill(rescaling, rescaling + rescaleWidth * 2, 0.0);

    isInvariant = new bool[numChar];
    for(int i = 0; i < numChar; i++){
        isInvariant[i] = rng.uniformRv() > 0.5;
    }

    activeT->updateAll();

    dpp->registerModel(this);
}

Model::~Model(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] activeCL;
    delete [] activeTP;
    delete [] isInvariant;
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

    invariantUpdate = false;
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

    invariantUpdate = false;
}

double Model::lnPrior(){
    return dpp->lnPrior() + tree->lnPrior() + rateMatrix->lnPrior();
}

double Model::updateInvariance() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    TreeObject* activeT = tree->getTree();
    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<int> assignments = dpp->getAssignments();

    for(int site = 0; site < numChar; site++){
        std::vector<double> likelihoods(2, 0.0);

        int invariant = (int)isInvariant[site];

        //First lets work with the invariant likelihood we have
        {
            int rIndex = activeT->getRoot()->getIndex();
            double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0) + site * stateSpace;
            std::vector<double> f = rateMatrix->getStationary();

            double lnL = 0.0;
            for(int i = 0; i < stateSpace; i++){
                lnL += pR[i]*f[i];
            }
            lnL = std::log(lnL);

            for(Node* n : poSeq){
                lnL += *(rescaling + (numChar * n->getIndex()) + site);
            }

            likelihoods[0] = lnL;
        }
        //Now lets calculate the likelihood for the other invariant state
        {
            int category = assignments[site];

            double* siteBuffer = new double[numNodes * stateSpace];
            std::fill(siteBuffer, siteBuffer + (numNodes * stateSpace), 1.0);

            double* rescaleBuffer = new double[numNodes];
            std::fill(rescaleBuffer, rescaleBuffer + numNodes, 0.0);

            for(Node* n : poSeq){
                int nIndex = n->getIndex();
                
                double* pNN = siteBuffer + (nIndex * stateSpace);
                if(n->getIsTip() == false){

                    std::set<Node*>& nNeighbors = n->getNeighbors();

                    for(Node* d : nNeighbors){
                        if(d != n->getAncestor()){
                            int dIndex = d->getIndex();
                            double* pN = pNN;
                            double* pD = siteBuffer + (dIndex * stateSpace);

                            Matrix<double> P = *(*transProb)(activeTP[dIndex], category*2 + (invariant ^ true), dIndex);
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

            likelihoods[1] = lnL;
        }

        double total = likelihoods[0] + likelihoods[1];
        double draw = rng.uniformRv() * total;

        if(draw > likelihoods[0])
            isInvariant[site] = !((bool)invariant);
    }

    invariantUpdate = true;
    return INFINITY;
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
        transProb->allocateQ(numCats*2);
        for(int i = 0; i < numCats; i++){
            double omega1 = categories[i].omega1;
            double omega2 = omega1 + categories[i].omega2;
            rateTaskflow.emplace([this, omega1, omega2, i](){
                for(int inv = 0; inv < 2; inv++)
                    transProb->updateQ(rateMatrix->Q(omega1, omega2, inv), i*2 + inv);
            });
        }

        executor.run(rateTaskflow).wait();
    }
    else if(dpp->isDirty()){
        activeT->updateAll();
        transProb->allocateQ(numCats*2);
        for(int i = 0; i < numCats; i++){
            if(categories[i].dirty){
                double omega1 = categories[i].omega1;
                double omega2 = omega1 + categories[i].omega2;
                for(int inv = 0; inv < 2; inv++)
                    transProb->updateQ(rateMatrix->Q(omega1, omega2, inv), i*2 + inv);
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
                    for(int i = 0; i < numCats; i++){
                        transProb->setProbs(activeTP[nIndex], i*2, nIndex, v);
                        transProb->setProbs(activeTP[nIndex], i*2 + 1, nIndex, v);
                    }
                }
                n->setNeedsTPUpdate(false);
            }
            if(n->getNeedsCLUpdate() == true || (invariantUpdate && n->getIsTip() == false)){
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
                        for(int cat = 0; cat < numCats; cat++){
                            transProbMatrices.push_back(*(*transProb)(activeTP[dIndex], cat*2, dIndex));
                            transProbMatrices.push_back(*(*transProb)(activeTP[dIndex], cat*2 + 1, dIndex));
                        }
                        for(int c = 0; c < numChar; c++){
                            int invariant = (int)isInvariant[c];
                            Matrix<double> P = transProbMatrices[assignments[c]*2 + invariant];
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
    double omega2 = omega1 + categories[category].omega2;
    
    for(int i = 0; i < 2; i++)
        transProb->updateQ(rateMatrix->Q(omega1, omega2, i), category*2 + i);

    for(Node* n : poSeq){
        int nIndex = n->getIndex();

        if(n != activeT->getRoot()){
            double v = activeT->getBranchLength(n);
            transProb->setProbs(activeTP[nIndex], category*2, nIndex, v);
            transProb->setProbs(activeTP[nIndex], category*2 + 1, nIndex, v);
        }
    }
}


double Model::testCategory(int site, int category, bool update){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*>&  poSeq = activeT->getPostOrderSeq();
    std::vector<Category> categories = dpp->getCategories();

    int invariant = (int)isInvariant[site];

    if(update) {
        double omega1 = categories[category].omega1;
        double omega2 = omega1 + categories[category].omega2;
        transProb->updateQ(rateMatrix->Q(omega1, omega2, invariant), category*2 + 1);
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
                transProb->setProbs(activeTP[nIndex], category*2 + invariant, nIndex, v);
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

                    Matrix<double> P = *(*transProb)(activeTP[dIndex], category*2 + invariant, dIndex);
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
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    std::string returnString = std::to_string(i);
    std::vector<Node*> tips = tree->getTree()->getTips();

    TreeObject* activeT = tree->getTree();
    std::vector<Node*>&  preOrderSeq = activeT->getPostOrderSeq();
    std::reverse(preOrderSeq.begin(), preOrderSeq.end());
    std::vector<int> assignments = dpp->getAssignments();
    std::vector<Category> categories = dpp->getCategories();
    int numCats = dpp->getNumCategories();

    int* reconstructedStates = new int[numNodes*numChar];
    double* reconstructedOmega = new double[numNodes*numChar];

    std::fill(reconstructedOmega, reconstructedOmega + numNodes*numChar, 0.0);

    int numJointDraws = 10;

    for(int d = 0; d < numJointDraws; d++){

        std::fill(reconstructedStates, reconstructedStates + numNodes*numChar, -1);

        tf::Taskflow phyloTaskflow;
        std::unordered_map<int, tf::Task> taskMap;

        Node* root = activeT->getRoot();
        taskMap.insert(std::make_pair(root->getIndex(), phyloTaskflow.emplace([this, root, &rng, reconstructedStates, reconstructedOmega, assignments, categories](){
            int rIndex = root->getIndex();
            double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
            int* reconstructedP = reconstructedStates + rIndex*numChar;
            double* omegaP = reconstructedOmega + rIndex*numChar;

            for(int c = 0; c < numChar; c++){
                double total = 0;
                for(int i = 0; i < stateSpace; i++){
                    total += pR[i];
                }

                double draw = rng.uniformRv() * total;

                double sum = 0;
                bool success = false;
                for(int i = 0; i < stateSpace; i++){
                    sum += pR[i];

                    if(sum >= draw){
                        *reconstructedP = i;
                        if(i < 61){
                            *omegaP += categories[assignments[c]].omega1;
                        }
                        else{
                            *omegaP += categories[assignments[c]].omega1 + categories[assignments[c]].omega2;
                        }
                        success = true;
                        break;
                    }
                }

                if(success == false)
                    Msg::error("Failed to reconstruct state!");

                pR += stateSpace;
                reconstructedP++;
                omegaP++;
            }
        })));

        for(Node* n : preOrderSeq){
            if(n != activeT->getRoot()){
                taskMap.insert(std::make_pair(n->getIndex(), phyloTaskflow.emplace([this, n, &rng, reconstructedStates, reconstructedOmega, categories, assignments](){
                    int nIndex = n->getIndex();
                    double* pN = (*postOrder)(nIndex, activeCL[nIndex], 0);

                    int* reconstructedP = reconstructedStates + nIndex*numChar;
                    double* omegaP = reconstructedOmega + nIndex*numChar;

                    int ancestorIndex = n->getAncestor()->getIndex();

                    for(int c = 0; c < numChar; c++){
                        int invariant = (int)isInvariant[c];
                        Matrix<double> P = *(*transProb)(activeTP[nIndex], assignments[c]*2 + invariant, nIndex);
                        int ancestorState = *(reconstructedStates + ancestorIndex*numChar + c);

                        double total = 0;
                        for(int i = 0; i < stateSpace; i++){
                            total += P(ancestorState, i) * pN[i];
                        }

                        double draw = rng.uniformRv() * total;

                        double sum = 0;
                        bool success = false;
                        for(int i = 0; i < stateSpace; i++){
                            sum += P(ancestorState, i) * pN[i];

                            if(sum >= draw){
                                *reconstructedP = i;
                                if(i < 61){
                                    *omegaP += categories[assignments[c]].omega1;
                                }
                                else{
                                    *omegaP += categories[assignments[c]].omega1 + categories[assignments[c]].omega2;
                                }
                                success = true;
                                break;
                            }
                        }

                        if(success == false)
                            Msg::error("Failed to reconstruct state!");

                        pN += stateSpace;
                        reconstructedP++;
                        omegaP++;
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

    for(Node* n : tips) {
        int index = n->getIndex();
        double* omegaP = reconstructedOmega + index*numChar;
        for(int c = 0; c < numChar; c++) {
            returnString += "\t" + std::to_string(*omegaP/numJointDraws);
            omegaP++;
        }
    }

    delete [] reconstructedStates;
    delete [] reconstructedOmega;

    return returnString + "\n";
}

std::string Model::invarHeader(){
    std::string returnString = "Iteration\tPosterior";
    for(int i = 0; i < numChar; i++)
        returnString += "\tInvariant[" + std::to_string(i) + "]";
    return returnString + "\n";
}

std::string Model::invarOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood);

    for(int i = 0; i < numChar; i++){
        returnString += "\t" + std::to_string(isInvariant[i]);
    }

    return returnString + "\n";
}