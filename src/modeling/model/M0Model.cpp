#include "M0Model.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/rate_matrices/M0Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <unordered_map>

#if TIME_PROFILE==1
#include <chrono>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

M0Model::M0Model(Settings* s, Alignment* a, TreeParameter* t, M0Matrix* m, tf::Executor& e) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), numChar(0), executor(e),
            branchLog(s->branchOutput), analysisLog(s->mcmcOutput), treeLog(s->treeOutput) {

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

    std::memcpy(rescaling + numNodes*numChar, rescaling, numNodes*numChar*sizeof(double));

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

    std::memcpy(rescaling, rescaling + numNodes*numChar, numNodes*numChar*sizeof(double));

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

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #endif

    if(rateMatrix->isDirty()){
        activeT->updateAll();
        transProb->updateQ(rateMatrix->Q(), 0);
    }

    
    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point rateTime = std::chrono::steady_clock::now();
    std::cout << "Rate computation was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(rateTime - begin).count() << "[milliseconds]" << std::endl;
    #endif

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

    
    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point probsTime = std::chrono::steady_clock::now();
    std::cout << "Probs computation was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(probsTime - rateTime).count() << "[milliseconds]" << std::endl;
    #endif

    for(Node* n : poSeq){
        if(n->getNeedsCLUpdate() == true){
            activeCL[n->getIndex()] ^= true; // Flip this ahead of time
        }
    }

    tf::Taskflow phyloTaskflow;
    
    int chunkSize = 25;
    for(int range = 0; range < (int)std::ceil((double)numChar / chunkSize); range++){
        int start = range * chunkSize;
        int end = start + chunkSize-1;
        end = std::min(end, numChar-1);

        phyloTaskflow.emplace([this, &poSeq, start, end](){

            // Avoid defining these within the tight inner loop and just treat them as working space
            #ifdef __AVX2__
            double tmp[4];
            __m256d pj;
            __m256d vj;
            #elif defined(__ARM_NEON__)
            float64x2_t pj;
            float64x2_t vj;
            float64x2_t prod;
            #endif

            int currentChunkSize = end - start + 1;
            for(Node* n : poSeq){
                int nIndex = n->getIndex();
                if(n->getNeedsCLUpdate() == true){
                    double* pNN = (*postOrder)(nIndex, activeCL[nIndex], 0) + start * stateSpace;
                    std::fill(pNN, pNN + (currentChunkSize * stateSpace), 1.0);

                    std::set<Node*>& nNeighbors = n->getNeighborRef();
                    for(Node* d : nNeighbors){
                        if(d != n->getAncestor()){
                            int dIndex = d->getIndex();
                            double* pN = pNN;
                            double* pD = (*postOrder)(dIndex, activeCL[dIndex], 0) + start * stateSpace;

                            const Matrix<double>& P = (*transProb)(activeTP[dIndex], 0, dIndex);
                            for(int c = 0; c < currentChunkSize; c++){
                                for(int i = 0; i < stateSpace; i++){
                                    double sum = 0.0;
                                    #ifdef __AVX2__
                                    int j = 0;

                                    // SIMD block - we will process multiples of 4 at a time
                                    __m256d sumVec = _mm256_setzero_pd();
                                    for (; j <= stateSpace - 4; j += 4) {
                                        pj = _mm256_loadu_pd(&P(i, j));
                                        vj = _mm256_loadu_pd(&pD[j]);
                                        sumVec = _mm256_fmadd_pd(pj, vj, sumVec); // += P(i,j) * pD(j)
                                    }

                                    _mm256_storeu_pd(tmp, sumVec); //Access all of the things we were doing simultaneous operations on and sum them
                                    sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];

                                    //In most of our models the state space is not perfectly divisible by 4
                                    for (; j < stateSpace; ++j) {
                                        sum += P(i, j) * pD[j];
                                    }
                                    #elif defined(__ARM_NEON__)
                                    int j = 0;
                                    
                                    // For the M series chips
                                    for (; j <= stateSpace - 2; j += 2) {
                                        pj = vld1q_f64(&P(i, j)); 
                                        vj = vld1q_f64(&pD[j]);
                                        prod = vmulq_f64(pj, vj);

                                        sum += vgetq_lane_f64(prod, 0) + vgetq_lane_f64(prod, 1);
                                    }

                                    // In some of our models the state space is not perfectly divisible by 2
                                    for (; j < stateSpace; ++j) {
                                        sum += P(i, j) * pD[j];
                                    }
                                    #else
                                    // In case our CPU really doesn't have any optimizations available.
                                    for(int j = 0; j < stateSpace; j++){
                                        sum += P(i, j) * pD[j];
                                    }
                                    #endif
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
                        pNN -= stateSpace;
                        for(int i = 0; i < stateSpace; i++){
                            *pNN /= max;
                            pNN++;
                        }
                        *rescalePointer = std::log(max);
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

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point pruneTime = std::chrono::steady_clock::now();
    std::cout << "Pruning was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(pruneTime - probsTime).count() << "[milliseconds]" << std::endl;
    #endif
}

void M0Model::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
}

std::vector<double> M0Model::getTunableParameterRecord() const {
    std::vector<double> record = {
        rateMatrix->getK(), rateMatrix->getOmega()
    };
    for(double entry : rateMatrix->getRawStationary())
        record.push_back(entry);
    for(double entry : tree->getTree()->getBranchLengths())
        record.push_back(entry);

    return record;
}

std::vector<double> M0Model::getTunableParameters() const {
    std::vector<double> returnVec(5, 0.0);
    returnVec[0] = tree->branchDelta;
    returnVec[1] = tree->treeAlpha;
    returnVec[2] = rateMatrix->stationaryAlpha;
    returnVec[3] = rateMatrix->kDelta;
    returnVec[4] = rateMatrix->omegaDelta;
    return returnVec;
}

void M0Model::setTunableParameters(const std::vector<double> & v){
    tree->branchDelta = v[0];
    tree->treeAlpha = v[1];
    rateMatrix->stationaryAlpha = v[2];
    rateMatrix->kDelta = v[3];
    rateMatrix->omegaDelta = v[4];
}

void M0Model::printAcceptanceRates(){
    std::cout << "Tree Acceptance Rate: " << tree->getTreeRate() << "\tBranch Acceptance Rate: " << tree->getBranchRate() <<
    "\tStationary Acceptance Rate: " << rateMatrix->getStationaryRate() << "\tK Acceptance Rate: " << rateMatrix->getKRate() <<
    "\tOmega Acceptance Rate: " << rateMatrix->getOmegaRate() << std::endl;
}

void M0Model::printTabular(int i){
    if(i == 0){
        std::string returnString = "Iteration\tPosterior\tLikelihood\tPrior\tK\tOmega";
        for(int i = 0; i < 61; i++)
            returnString += "\tPi[" + std::to_string(i) + "]";

        std::cout << returnString << std::endl;
    }
    else{
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                                std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" +
                                std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getOmega());
        std::vector<double> stationary = rateMatrix->getRawStationary();
        for(double i : stationary){
            returnString += "\t" + std::to_string(i);
        }

        std::cout << returnString << std::endl;
    }
}

void M0Model::writeLogHeaders(){
    if(analysisLog != ""){
        std::string tabHeader = "Iteration\tPosterior\tLikelihood\tPrior\tK\tOmega";
        for(int i = 0; i < 61; i++)
            tabHeader += "\tPi[" + std::to_string(i) + "]";
        tabHeader += "\n";
        std::ofstream outFile(analysisLog, std::ios::out);
        outFile << tabHeader;
    }

    if(treeLog != ""){
        std::string treeHeader =  "Iteration\tPosterior\tTree\n";
        std::ofstream outFile(treeLog, std::ios::out);
        outFile << treeHeader;
    }

    if(branchLog != ""){
        std::string branchHeader = "Iteration\tPosterior";
        TreeObject* treeObj = tree->getTree();
        std::vector<Node*> poSeq = treeObj->getPostOrderSeq();
        std::vector<Node*> nodes;
        for(Node* n : poSeq)
            nodes.push_back(n);
        std::sort(nodes.begin(), nodes.end(), [](const Node* a, const Node* b) {
            return a->getIndex() > b->getIndex();
        });

        for(Node* n : nodes){
            if(n != treeObj->getRoot()){
                branchHeader += "\tBranch[" + std::to_string(n->getIndex()) + "]";
            }
        }
        branchHeader += "\n";

        std::ofstream outFile(branchLog, std::ios::out);
        outFile << branchHeader;
    }
}

void M0Model::writeLogData(int i) {
    if(analysisLog != ""){
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                                std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" +
                                std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getOmega());
        std::vector<double> stationary = rateMatrix->getRawStationary();
        for(double i : stationary){
            returnString += "\t" + std::to_string(i);
        }
        returnString += "\n";

        std::ofstream outFile(analysisLog, std::ios::app);
        outFile << returnString;
    }

    if(treeLog != ""){
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" + tree->writeNewick() + "\n";
        
        std::ofstream outFile(treeLog, std::ios::app);
        outFile << returnString;
    }

    if(branchLog != ""){
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood);
        TreeObject* treeObj = tree->getTree();
        std::vector<Node*> poSeq = treeObj->getPostOrderSeq();
        std::vector<Node*> nodes;
        for(Node* n : poSeq)
            nodes.push_back(n);
        std::sort(nodes.begin(), nodes.end(), [](const Node* a, const Node* b) {
            return a->getIndex() > b->getIndex();
        });

        for(Node* n : nodes){
            if(n != treeObj->getRoot()){
                returnString += "\t" + std::to_string(treeObj->getBranchLength(n));
            }
        }

        returnString =+ "\n";

        std::ofstream outFile(branchLog, std::ios::app);
        outFile << returnString;
    }
}