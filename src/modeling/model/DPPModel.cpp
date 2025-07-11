#include "DPPModel.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/DPPMatrix.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
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

DPPModel::DPPModel(Settings s, Alignment* a, TreeParameter* t, DPPMatrix* m, DirichletProcessPrior* d) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0),
            dpp(d), numChar(0) {

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
    postOrder = new ConditionalLikelihood(aln, numNodes, 1, stateSpace);
    transProb = new TransitionProbability(numNodes, stateSpace);

    int rescaleWidth = numNodes*numChar;
    rescaling = new double[rescaleWidth * 2];
    std::fill(rescaling, rescaling + rescaleWidth * 2, 0.0);

    activeT->updateAll();

    dpp->registerModel(this);
}

DPPModel::~DPPModel(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] activeCL;
    delete [] activeTP;
}

void DPPModel::accept() {
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
    if(dpp->isDirty()){
        dpp->accept();
        dpp->clean();
    }

    transProb->accept();
}

void DPPModel::reject() {
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
    if(dpp->isDirty()){
        dpp->reject();
        dpp->clean();
    }

    transProb->reject();
}

double DPPModel::lnPrior(){
    return dpp->lnPrior() + tree->lnPrior() + rateMatrix->lnPrior();
}

void DPPModel::regenerateLikelihood(){
    TreeObject* activeT = tree->getTree();

    const std::vector<Node*> poSeq = activeT->getPostOrderSeq();
    std::vector<int> assignments = dpp->getAssignments();
    std::vector<Category> categories = dpp->getCategories();
    int numCats = dpp->getNumCategories();

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #endif

    tf::Taskflow rateTaskflow;

    if(rateMatrix->isDirty()){
        activeT->updateAll();
        transProb->allocateQ(numCats);
        for(int i = 0; i < numCats; i++){
            double omega1 = categories[i].omega1;
            double omega2 = categories[i].omega2;
            rateTaskflow.emplace([this, i, omega1, omega2](){
                transProb->updateQ(rateMatrix->Q(omega1, omega2), i);
            });
        }
    }
    else if(dpp->isDirty()){ // Only spend time updating the dirty ones
        activeT->updateAll();
        transProb->allocateQ(numCats);
        for(int i = 0; i < numCats; i++){
            if(categories[i].dirty){
                double omega1 = categories[i].omega1;
                double omega2 = categories[i].omega2;
                rateTaskflow.emplace([this, i, omega1, omega2](){
                    transProb->updateQ(rateMatrix->Q(omega1, omega2), i);
                });
            }
        }
    }

    executor.run(rateTaskflow).wait();

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
                probsTaskflow.emplace([this, numCats, nIndex, v, activeIndex](){
                    for(int i = 0; i < numCats; i++){
                        transProb->setProbs(activeIndex, i, nIndex, v);
                    }
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
    
    int chunkSize = 50;
    for(int range = 0; range < (int)std::ceil((double)numChar / chunkSize); range++){
        int start = range * chunkSize;
        int end = start + chunkSize-1;
        end = std::min(end, numChar-1);

        phyloTaskflow.emplace([this, &poSeq, numCats, &assignments, start, end](){
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

                            for(int c = 0; c < currentChunkSize; c++){
                                const Matrix<double>& P = (*transProb)(activeTP[dIndex], assignments[c + start], dIndex);
                                for(int i = 0; i < stateSpace; i++){
                                    double sum = 0.0;
                                    #ifdef __AVX2__
                                    int j = 0;

                                    // SIMD block - we will process multiples of 4 at a time
                                    __m256d sumVec = _mm256_setzero_pd();
                                    for (; j <= stateSpace - 4; j += 4) {
                                        __m256d pj = _mm256_loadu_pd(&P(i, j));
                                        __m256d vj = _mm256_loadu_pd(&pD[j]);
                                        sumVec = _mm256_fmadd_pd(pj, vj, sumVec); // += P(i,j) * pD(j)
                                    }

                                    double tmp[4];
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
                                        float64x2_t pj = vld1q_f64(&P(i, j)); 
                                        float64x2_t vj = vld1q_f64(&pD[j]);
                                        float64x2_t prod = vmulq_f64(pj, vj);

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

    #if LOGGING == 1
    std::cout << "Non-rescaled likelihood: " << lnL << std::endl;
    #endif

    double* rescalePointer = rescaling;
    for(int i = 0, len = numNodes * numChar; i < len; i++){
        lnL += *rescalePointer;
        rescalePointer++;
    }

    #if LOGGING == 1
    std::cout << "Rescaled likelihood: " << lnL << std::endl;
    #endif

    currentLikelihood = lnL;

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point pruneTime = std::chrono::steady_clock::now();
    std::cout << "Pruning was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(pruneTime - probsTime).count() << "[milliseconds]" << std::endl;
    #endif
}

void DPPModel::regenerateTransitionProbs(int site, int category){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*> poSeq = activeT->getPostOrderSeq();
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


double DPPModel::testCategory(int site, int category, bool update){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*> poSeq = activeT->getPostOrderSeq();
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

    // Some setup
    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        if(update){
            if(n != activeT->getRoot()){
                double v = activeT->getBranchLength(n);
                transProb->setProbs(activeTP[nIndex], category, nIndex, v);
            }
        }

        if(n->getIsTip()){
            double* dataPointer = (*postOrder)(nIndex, activeCL[nIndex], 0) + site*stateSpace;
            double* sitePointer = siteBuffer + nIndex*stateSpace;
            for(int i = 0; i < stateSpace; i++){
                *sitePointer = *dataPointer;
                dataPointer++;
                sitePointer++;
            }
        }
    }

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

                    Matrix<double> P = (*transProb)(activeTP[dIndex], category, dIndex);
                    for(int i = 0; i < stateSpace; i++){
                        double sum = 0.0;
                        #ifdef __AVX2__
                        int j = 0;

                        // SIMD block - we will process multiples of 4 at a time
                        __m256d sumVec = _mm256_setzero_pd();
                        for (; j <= stateSpace - 4; j += 4) {
                            __m256d pj = _mm256_loadu_pd(&P(i, j));
                            __m256d vj = _mm256_loadu_pd(&pD[j]);
                            sumVec = _mm256_fmadd_pd(pj, vj, sumVec); // += P(i,j) * pD(j)
                        }

                        double tmp[4];
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
                            float64x2_t pj = vld1q_f64(&P(i, j)); 
                            float64x2_t vj = vld1q_f64(&pD[j]);
                            float64x2_t prod = vmulq_f64(pj, vj);

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
                for(int i = 0; i < stateSpace; i++){
                    *pNN /= max;
                    pNN++;
                }
                *rescalePointer = std::log(max);
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

void DPPModel::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
    dpp->tune();
}

std::string DPPModel::tabularHeader(){
    std::string returnString = "Iteration\tPosterior\tLikelihood\tPrior\tK\tR";
    for(int i = 0; i < 61; i++){
        returnString += "\tPi[" + std::to_string(i) + "]";
    }

    return returnString + "\n";
}

std::string DPPModel::tabularOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                               std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" +
                               std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getR());
    std::vector<double> stationary = rateMatrix->getRawStationary();
    for(double i : stationary){
        returnString += "\t" + std::to_string(i);
    }

    return returnString + "\n";
}

std::string DPPModel::treeHeader(){
    return "Iteration\tPosterior\tTree\n";
}

std::string DPPModel::treeOut(int i){
    return std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" + tree->writeNewick() + "\n";
}

std::string DPPModel::dppHeader(){
    std::string returnString = "Iteration\tPosterior\tCategoryCount";
    for(int i = 0; i < numChar; i++)
        returnString += "\tOmega1[" + std::to_string(i) + "]" + "\tOmegaIncrement[" + std::to_string(i) + "]";
    return returnString + "\n";
}

std::string DPPModel::dppOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t";
    std::vector<Category> categories = dpp->getCategories();
    returnString += std::to_string(categories.size());
    std::vector<int> assignments = dpp->getAssignments();
    for(int c : assignments){
        returnString += "\t" + std::to_string(categories[c].omega1) + "\t" + std::to_string(categories[c].omegaIncrement);
    }


    return returnString + "\n";
}

std::string DPPModel::tipsHeader(){
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

std::string DPPModel::ancestralHeader(){
    std::string returnString = "Iteration";

    for(int i = 0; i < numNodes; i++) {
        for(int c = 0; c < numChar; c++){
            returnString += "\t" + std::to_string(i) + "[" + std::to_string(c) + "]";
        }
    }

    return returnString + "\n";
}

std::tuple<std::string, std::string> DPPModel::reconstructionOut(int i){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    std::string tipString = std::to_string(i);
    std::string ancestralString = std::to_string(i);
    std::vector<Node*> tips = tree->getTree()->getTips();

    TreeObject* activeT = tree->getTree();
    std::vector<Node*> preOrderSeq = activeT->getPostOrderSeq();
    std::reverse(preOrderSeq.begin(), preOrderSeq.end());
    std::vector<int> assignments = dpp->getAssignments();
    std::vector<Category> categories = dpp->getCategories();
    int numCats = dpp->getNumCategories();

    std::vector<double> dNdS1;
    std::vector<double> dNdS2;
    for(Category c : categories){
        auto ratio_pair = rateMatrix->dNdS(c.omega1, c.omega2);
        dNdS1.push_back(ratio_pair.first);
        dNdS2.push_back(ratio_pair.second);
    }

    int* reconstructedStates = new int[numNodes*numChar];
    double* reconstructeddNdS = new double[numNodes*numChar];

    std::fill(reconstructeddNdS, reconstructeddNdS + numNodes*numChar, 0.0);

    int numJointDraws = 10;

    for(int d = 0; d < numJointDraws; d++){

        std::fill(reconstructedStates, reconstructedStates + numNodes*numChar, -1);

        tf::Taskflow phyloTaskflow;
        std::unordered_map<int, tf::Task> taskMap;

        Node* root = activeT->getRoot();
        taskMap.insert(std::make_pair(root->getIndex(), phyloTaskflow.emplace([this, root, &rng, reconstructedStates, reconstructeddNdS, assignments, dNdS1, dNdS2](){
            int rIndex = root->getIndex();
            double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
            int* reconstructedP = reconstructedStates + rIndex*numChar;
            double* dNdSP = reconstructeddNdS + rIndex*numChar;

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
                            *dNdSP += dNdS1[assignments[c]];
                        }
                        else{
                            *dNdSP += dNdS2[assignments[c]];
                        }
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
                taskMap.insert(std::make_pair(n->getIndex(), phyloTaskflow.emplace([this, n, &rng, reconstructedStates, reconstructeddNdS, assignments, dNdS1, dNdS2](){
                    int nIndex = n->getIndex();
                    double* pN = (*postOrder)(nIndex, activeCL[nIndex], 0);

                    int* reconstructedP = reconstructedStates + nIndex*numChar;
                    double* dNdSP = reconstructeddNdS + nIndex*numChar;

                    int ancestorIndex = n->getAncestor()->getIndex();

                    for(int c = 0; c < numChar; c++){
                        Matrix<double> P = (*transProb)(activeTP[nIndex], assignments[c], nIndex);
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
                                    *dNdSP += dNdS1[assignments[c]];
                                }
                                else{
                                    *dNdSP += dNdS2[assignments[c]];
                                }
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
        double* dNdSP = reconstructeddNdS + index*numChar;
        for(int c = 0; c < numChar; c++) {
            tipString += "\t" + std::to_string(*dNdSP/numJointDraws);
            dNdSP++;
        }
    }

    for(int i = 0; i < numNodes; i++) {
        double* dNdSP = reconstructeddNdS + i*numChar;
        for(int c = 0; c < numChar; c++) {
            ancestralString += "\t" + std::to_string(*dNdSP/numJointDraws);
            dNdSP++;
        }
    }

    delete [] reconstructedStates;
    delete [] reconstructeddNdS;

    return std::make_tuple(tipString + "\n", ancestralString + "\n");
}

std::string DPPModel::branchHeader(){
    std::string returnString = "Iteration\tPosterior";
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
            returnString += "\tBranch[" + std::to_string(n->getIndex()) + "]";
        }
    }
    return returnString + "\n";
}

std::string DPPModel::branchOut(int i){
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

    return returnString + "\n";
}
