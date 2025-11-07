#include "RJModel.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/rate_matrices/RJMatrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <unordered_map>
#include <chrono>
#include <fstream>

#if TIME_PROFILE==1
#include <chrono>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

RJModel::RJModel(Settings* s, Alignment* a, TreeParameter* t, RJMatrix* m, tf::Executor& e) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), numChar(0), executor(e),
            branchLog(s->branchOutput), tipsLog(s->tipsOutput), ancestralLog(s->ancestralStatesOutput), treeLog(s->treeOutput),
            analysisLog(s->mcmcOutput)  {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    TreeObject* activeT = tree->getTree();
    stateSpace = 305;
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

RJModel::~RJModel(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] activeCL;
    delete [] activeTP;
}

void RJModel::accept() {
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

void RJModel::reject() {
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

double RJModel::lnPrior(){
    double prior = tree->lnPrior() + rateMatrix->lnPrior();
    int activeOmegas = rateMatrix->getActiveOmegas();

    // Prior over the number of active omegas
    if(activeOmegas == 1){
        prior += std::log(0.25);
    }
    else if(activeOmegas == 2){
        prior += std::log(0.25);
    }
    else if(activeOmegas == 3){
        prior += std::log(0.20);
    }
    else if(activeOmegas == 4){
        prior += std::log(0.15);
    }
    else{
        prior += std::log(0.15);
    }


    return prior;
}

void RJModel::regenerateLikelihood(){
    TreeObject* activeT = tree->getTree();

    const std::vector<Node*> poSeq = activeT->getPostOrderSeq();

    int numClasses = rateMatrix->getActiveOmegas();
    int activeSubspace = numClasses * 61;

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #endif

    if(rateMatrix->isDirty()){
        activeT->updateAll();
        transProb->updateQ(rateMatrix->Q(), 0); // This should automatically adjust to the size of the matrix
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
                probsTaskflow.emplace([this, nIndex, v, activeIndex, numClasses, activeSubspace](){
                    transProb->setProbs(activeIndex, 0, nIndex, activeSubspace, v);
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

        phyloTaskflow.emplace([this, &poSeq, start, end, numClasses, activeSubspace](){

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
                                for(int i = 0; i < activeSubspace; i++){
                                    double sum = 0.0;
                                    #ifdef __AVX2__
                                    int j = 0;

                                    // SIMD block - we will process multiples of 4 at a time
                                    __m256d sumVec = _mm256_setzero_pd();
                                    for (; j <= activeSubspace - 4; j += 4) {
                                        pj = _mm256_loadu_pd(&P(i, j));
                                        vj = _mm256_loadu_pd(&pD[j]);
                                        sumVec = _mm256_fmadd_pd(pj, vj, sumVec); // += P(i,j) * pD(j)
                                    }

                                    _mm256_storeu_pd(tmp, sumVec); //Access all of the things we were doing simultaneous operations on and sum them
                                    sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];

                                    //In most of our models the state space is not perfectly divisible by 4
                                    for (; j < activeSubspace; ++j) {
                                        sum += P(i, j) * pD[j];
                                    }
                                    #elif defined(__ARM_NEON__)
                                    int j = 0;
                                    
                                    // For the M series chips
                                    for (; j <= activeSubspace - 2; j += 2) {
                                        pj = vld1q_f64(&P(i, j)); 
                                        vj = vld1q_f64(&pD[j]);
                                        prod = vmulq_f64(pj, vj);

                                        sum += vgetq_lane_f64(prod, 0) + vgetq_lane_f64(prod, 1);
                                    }

                                    // In some of our models the state space is not perfectly divisible by 2
                                    for (; j < activeSubspace; ++j) {
                                        sum += P(i, j) * pD[j];
                                    }
                                    #else
                                    // In case our CPU really doesn't have any optimizations available.
                                    for(int j = 0; j < activeSubspace; j++){
                                        sum += P(i, j) * pD[j];
                                    }
                                    #endif
                                    (*pN) *= sum;
                                    pN++;
                                }
                                pN+=(5-numClasses)*61; // Shift unused classes
                                pD+=stateSpace;
                            }
                        }
                    }
                    double* rescalePointer = rescaling + (numChar * nIndex) + start;
                    std::fill(rescalePointer, rescalePointer + currentChunkSize, 0.0);

                    for(int c = 0; c < currentChunkSize; c++){
                        double max = *pNN;
                        pNN++;
                        for(int i = 1; i < activeSubspace; i++){
                            if(*pNN > max)
                                max = *pNN;
                            pNN++;
                        }
                        pNN -= activeSubspace;
                        for(int i = 0; i < activeSubspace; i++){
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
        for(int i = 0; i < numClasses * 61; i++){
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

void RJModel::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
}

std::vector<double> RJModel::getTunableParameterRecord() const {
    std::vector<double> record = {
        rateMatrix->getK(), rateMatrix->getOmega(0), rateMatrix->getOmega(1), 
        rateMatrix->getOmega(2), rateMatrix->getOmega(3), rateMatrix->getOmega(4),
        rateMatrix->getR()
    };
    for(double entry : rateMatrix->getRawStationary())
        record.push_back(entry);
    for(double entry : tree->getTree()->getBranchProportions())
        record.push_back(entry);

    return record;
}

std::vector<double> RJModel::getTunableParameters() const {
    std::vector<double> returnVec(6, 0.0);
    returnVec[0] = tree->branchAlpha;
    returnVec[1] = tree->treeDelta;
    returnVec[2] = rateMatrix->stationaryAlpha;
    returnVec[3] = rateMatrix->kDelta;
    returnVec[4] = rateMatrix->omegaDelta;
    returnVec[5] = rateMatrix->rDelta;
    return returnVec;
}

void RJModel::setTunableParameters(const std::vector<double> & v){
    tree->branchAlpha = v[0];
    tree->treeDelta = v[1];
    rateMatrix->stationaryAlpha = v[2];
    rateMatrix->kDelta = v[3];
    rateMatrix->omegaDelta = v[4];
    rateMatrix->rDelta = v[5];
}

void RJModel::printAcceptanceRates() {
    std::cout << "Tree Acceptance Rate: " << tree->getTreeRate() << "\tBranch Acceptance Rate: " << tree->getBranchRate() <<
    "\tStationary Acceptance Rate: " << rateMatrix->getStationaryRate() << "\tK Acceptance Rate: " << rateMatrix->getKRate() <<
    "\tOmega Acceptance Rate: " << rateMatrix->getOmegaRate() << "\tR Acceptance Rate: " << rateMatrix->getRRate() << std::endl;
}

void RJModel::printTabular(int i) {
    if(i == 0){
        std::string returnString = "Iteration\tPosterior\tLikelihood\tPrior\tTreeLength\tOmegaCount\tK\tOmega\tOmegaIncrement1\tOmegaIncrement2\tOmegaIncrement3\tOmegaIncrement4\tR";
        for(int i = 0; i < 61; i++)
            returnString += "\tPi[" + std::to_string(i) + "]";
            

        std::cout << returnString << std::endl;
    }
    else{
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                                std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" + std::to_string(tree->getTree()->getTreeLength()) + "\t" +
                                std::to_string(rateMatrix->getActiveOmegas()) + "\t" +
                                std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getOmega(0)) + "\t" +
                                std::to_string(rateMatrix->getOmega(1)) + "\t" + std::to_string(rateMatrix->getOmega(2)) + "\t" +
                                std::to_string(rateMatrix->getOmega(3)) + "\t" + std::to_string(rateMatrix->getOmega(4)) + "\t" +
                                std::to_string(rateMatrix->getR());
        std::vector<double> stationary = rateMatrix->getRawStationary();
        for(double i : stationary){
            returnString += "\t" + std::to_string(i);
        }

        std::cout << returnString << std::endl;
    }
}

void RJModel::writeLogHeaders() {
    if(analysisLog != ""){
        std::string tabHeader = "Iteration\tPosterior\tLikelihood\tPrior\tTreeLength\tOmegaCount\tK\tOmega\tOmegaIncrement1\tOmegaIncrement2\tOmegaIncrement3\tOmegaIncrement4\tR";
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

    if(tipsLog != ""){
        std::string tipHeader = "Iteration";
        std::vector<Node*> tips = tree->getTree()->getTips();
        for(Node* n : tips){
            std::string name = n->getName();
            for(int c = 0; c < numChar; c++){
                tipHeader += "\t" + name + "[" + std::to_string(c) + "]";
            }
        }
        tipHeader += "\n";

        std::ofstream outFile(tipsLog, std::ios::out);
        outFile << tipHeader;
    }

    if(ancestralLog != ""){
        std::string ancestralHeader = "Iteration";
        for(int i = 0; i < numNodes; i++) {
            for(int c = 0; c < numChar; c++){
                ancestralHeader += "\t" + std::to_string(i) + "[" + std::to_string(c) + "]";
            }
        }
        ancestralHeader += "\n";

        std::ofstream outFile(ancestralLog, std::ios::out);
        outFile << ancestralHeader;
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

void RJModel::writeLogData(int i) {
    if(analysisLog != ""){
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                                std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" + std::to_string(tree->getTree()->getTreeLength()) + "\t" +
                                std::to_string(rateMatrix->getActiveOmegas()) + "\t" +
                                std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getOmega(0)) + "\t" +
                                std::to_string(rateMatrix->getOmega(1)) + "\t" + std::to_string(rateMatrix->getOmega(2)) + "\t" +
                                std::to_string(rateMatrix->getOmega(3)) + "\t" + std::to_string(rateMatrix->getOmega(4)) + "\t" +
                                std::to_string(rateMatrix->getR());
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

    if(tipsLog != "" || ancestralLog != ""){
        RandomVariable& rng = RandomVariable::randomVariableInstance();

        std::vector<Node*> tips = tree->getTree()->getTips();

        TreeObject* activeT = tree->getTree();
        std::vector<Node*> preOrderSeq = activeT->getPostOrderSeq();
        std::reverse(preOrderSeq.begin(), preOrderSeq.end());

        std::array<double, 5> dNdS = rateMatrix->dNdS();
        int numClasses = rateMatrix->getActiveOmegas();
        int activeSubspace = 61*numClasses;

        int* reconstructedStates = new int[numNodes*numChar];
        double* reconstructeddNdS = new double[numNodes*numChar];

        std::fill(reconstructeddNdS, reconstructeddNdS + numNodes*numChar, 0.0);
        std::fill(reconstructedStates, reconstructedStates + numNodes*numChar, -1);

        tf::Taskflow phyloTaskflow;
        std::unordered_map<int, tf::Task> taskMap;

        Node* root = activeT->getRoot();
        taskMap.insert(std::make_pair(root->getIndex(), phyloTaskflow.emplace([this, root, &rng, reconstructedStates, reconstructeddNdS, dNdS, numClasses, activeSubspace](){
            int rIndex = root->getIndex();
            double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
            int* reconstructedP = reconstructedStates + rIndex*numChar;
            double* dNdSP = reconstructeddNdS + rIndex*numChar;

            for(int c = 0; c < numChar; c++){
                double total = 0;
                for(int i = 0; i < activeSubspace; i++){
                    total += pR[i];
                }

                double draw = rng.uniformRv() * total;

                double sum = 0;
                bool success = false;
                for(int i = 0; i < activeSubspace; i++){
                    sum += pR[i];

                    if(sum >= draw){
                        *reconstructedP = i;
                        *dNdSP = dNdS[(int)(i/61.0)];
                        success = true;
                        break;
                    }
                }

                if(success == false)
                    Msg::error("Failed to reconstruct state!");

                pR += stateSpace;
                reconstructedP++;
                dNdSP++;
            }
        })));

        for(Node* n : preOrderSeq){
            if(n != activeT->getRoot()){
                taskMap.insert(std::make_pair(n->getIndex(), phyloTaskflow.emplace([this, n, &rng, reconstructedStates, reconstructeddNdS, dNdS, numClasses, activeSubspace](){
                    int nIndex = n->getIndex();
                    double* pN = (*postOrder)(nIndex, activeCL[nIndex], 0);

                    int* reconstructedP = reconstructedStates + nIndex*numChar;
                    double* dNdSP = reconstructeddNdS + nIndex*numChar;

                    int ancestorIndex = n->getAncestor()->getIndex();

                    for(int c = 0; c < numChar; c++){
                        Matrix<double> P = (*transProb)(activeTP[nIndex], 0, nIndex);
                        int ancestorState = *(reconstructedStates + ancestorIndex*numChar + c);

                        double total = 0;
                        for(int i = 0; i < activeSubspace; i++){
                            total += P(ancestorState, i) * pN[i];
                        }

                        double draw = rng.uniformRv() * total;

                        double sum = 0;
                        bool success = false;
                        for(int i = 0; i < activeSubspace; i++){
                            sum += P(ancestorState, i) * pN[i];

                            if(sum >= draw){
                                *reconstructedP = i;
                                *dNdSP = dNdS[(int)(i/61.0)];
                                success = true;
                                break;
                            }
                        }

                        if(success == false)
                            Msg::error("Failed to reconstruct state!");

                        pN += stateSpace;
                        reconstructedP++;
                        dNdSP++;
                    }
                })));
            }
        }

        for(Node* n : preOrderSeq){
            if(n->getIsTip() == false){
                for(Node* d : n->getNeighborRef()){
                    if(d != n->getAncestor()){
                        taskMap.at(n->getIndex()).precede(taskMap.at(d->getIndex()));
                    }
                }
            }
        }

        executor.run(phyloTaskflow).wait();

        if(tipsLog != ""){
            std::string tipString = std::to_string(i);
            for(Node* n : tips) {
                int index = n->getIndex();
                double* dNdSP = reconstructeddNdS + index*numChar;
                for(int c = 0; c < numChar; c++) {
                    tipString += "\t" + std::to_string(*dNdSP);
                    dNdSP++;
                }
            }

            std::ofstream outFile(tipsLog, std::ios::app);
            outFile << tipString;
        }

        if(ancestralLog != ""){
            std::string ancestralString = std::to_string(i);
            for(int i = 0; i < numNodes; i++) {
                double* dNdSP = reconstructeddNdS + i*numChar;
                for(int c = 0; c < numChar; c++) {
                    ancestralString += "\t" + std::to_string(*dNdSP);
                    dNdSP++;
                }
            }

            std::ofstream outFile(ancestralLog, std::ios::app);
            outFile << ancestralString;
        }

        delete [] reconstructedStates;
        delete [] reconstructeddNdS;
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