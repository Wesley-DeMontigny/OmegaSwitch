#include "DPCMMModel.hpp"
#include "misc/RandomVariable.hpp"
#include "misc/Alignment.hpp"
#include "misc/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/rate_matrices/DPCMMMatrix.hpp"
#include "misc/RandomVariable.hpp"
#include "misc/Settings.hpp"
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

DPCMMModel::DPCMMModel(Settings* s, Alignment* a, TreeParameter* t, DPCMMMatrix* m, tf::Executor& e) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), numChar(0), executor(e), numGibbsUpdate(s->numGibbs),
            tipsLog(s->tipsOutput), ancestralLog(s->ancestralStatesOutput), treeLog(s->treeOutput),
            analysisLog(s->mcmcOutput), dppLog(s->dppOutput), omegaLambda(s->omegaLambda), samplingPower(1.0) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    TreeObject* activeT = tree->getTree();
    stateSpace = 183;
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

DPCMMModel::~DPCMMModel(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] activeCL;
    delete [] activeTP;
}

void DPCMMModel::accept() {
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
        transProb->accept(); // TransProb only needs to reject/accept when the eigenvalues have changed
    }

    transProb->accept();
}

void DPCMMModel::reject() {
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
        transProb->reject(); // TransProb only needs to reject/accept when the eigenvalues have changed
    }
}

double DPCMMModel::lnPrior(){
    double prior = tree->lnPrior() + rateMatrix->lnPrior();

    return prior;
}

void DPCMMModel::regenerateLikelihood(){
    #ifdef SAMPLE_PRIOR
        return;
    #endif

    TreeObject* activeT = tree->getTree();

    const std::vector<Node*> poSeq = activeT->getPostOrderSeq();
    std::vector<int> assignments = rateMatrix->getAssignments();
    std::vector<Category> categories = rateMatrix->getCategories();
    int numCats = rateMatrix->getNumCategories();
    int numClasses = rateMatrix->getActiveOmegas();
    int activeSubspace = 61 * numClasses;

    #if TIME_PROFILE==1
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #endif

    tf::Taskflow rateTaskflow;

    if(rateMatrix->isDirty()){
        activeT->updateAll();
        transProb->allocateQ(numCats);
        for(int i = 0; i < numCats; i++){
            rateTaskflow.emplace([this, i](){
                transProb->updateQ(rateMatrix->Q(i), i);
            });
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
                probsTaskflow.emplace([this, numCats, nIndex, v, activeIndex, activeSubspace](){
                    for(int i = 0; i < numCats; i++){
                        transProb->setProbs(activeIndex, i, nIndex, activeSubspace, v);
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
    
    int chunkSize = (int)std::ceil((double)numChar/(double)executor.num_workers());
    for(int range = 0; range < (int)std::ceil((double)numChar / (double)chunkSize); range++){
        int start = range * chunkSize;
        int end = std::min(start + chunkSize, numChar) - 1;

        phyloTaskflow.emplace([this, &poSeq, numCats, &assignments, numClasses, activeSubspace, start, end](){

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

                            for(int c = 0; c < currentChunkSize; c++){
                                const Matrix<double>& P = (*transProb)(activeTP[dIndex], assignments[c + start], dIndex);
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
                                pN+=(3-numClasses)*61;
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
        for(int i = 0; i < activeSubspace; i++){
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

void DPCMMModel::regenerateTransitionProbs(int category){
    TreeObject* activeT = tree->getTree();

    std::vector<Node*> poSeq = activeT->getPostOrderSeq();
    int numClasses = rateMatrix->getActiveOmegas();
    int activeSubspace = numClasses * 61;
    
    transProb->updateQ(rateMatrix->Q(category), category);

    for(Node* n : poSeq){
        int nIndex = n->getIndex();

        if(n != activeT->getRoot()){
            double v = activeT->getBranchLength(n);
            transProb->setProbs(activeTP[nIndex], category, nIndex, activeSubspace, v);
        }
    }
}

double DPCMMModel::updateDPP(){
    if(rateMatrix->isAssignmentFixed()){
        return -1 * INFINITY;
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    TreeObject* activeT = tree->getTree();
    std::vector<Node*> poSeq = activeT->getPostOrderSeq();

    rateMatrix->dirty();

    int numClasses = rateMatrix->getActiveOmegas();
    int activeSubspace = numClasses * 61;
    int numAux = 5;
    int bufferSize = ((2 * numAux) + rateMatrix->getNumCategories());

    double alphaSplit = std::log(rateMatrix->getAlpha()/(double)numAux);

    double* tempCLBuffer = new double[numNodes * bufferSize * stateSpace];
    double* tempRescaleBuffer = new double[numNodes * bufferSize];
    std::fill(tempCLBuffer, tempCLBuffer + numNodes * bufferSize * stateSpace, 0.0);
    std::fill(tempRescaleBuffer, tempRescaleBuffer + numNodes * bufferSize, 0.0);

    #ifdef __AVX2__
    double tmp[4];
    __m256d pj;
    __m256d vj;
    #elif defined(__ARM_NEON__)
    float64x2_t pj;
    float64x2_t vj;
    float64x2_t prod;
    #endif
    
    for(int iter = 0; iter < numGibbsUpdate; iter++) {
        int randomSite = (int)(rng.uniformRv() * numChar);

        auto deleted = rateMatrix->unassignMember(randomSite); // This will also return an integer if it needs to be deleted
        if(deleted.has_value()){
           transProb->deleteQ(deleted.value());
        }

        int numCurrentCats = rateMatrix->getNumCategories();
        int numTestableCats = numCurrentCats + numAux;
        std::vector<Category> currentCategories = rateMatrix->getCategories();
        std::vector<std::array<double, 3>> newOmegas;

        for(int n = 0; n < numAux; n++){
            newOmegas.emplace_back(std::array<double,3>{
                Probability::Exponential::rv(&rng, omegaLambda),
                Probability::Exponential::rv(&rng, omegaLambda),
                Probability::Exponential::rv(&rng, omegaLambda)
            });
        }

        #ifndef SAMPLE_PRIOR

        transProb->allocateQ(numCurrentCats + (numAux * 2));

        tf::Taskflow transProbTasks;

        for(int c = numCurrentCats; c < numTestableCats; c++){
            transProbTasks.emplace([this, &poSeq, c, numCurrentCats, activeSubspace, &newOmegas, activeT](){
                transProb->updateQ(rateMatrix->Q(newOmegas[c - numCurrentCats]), c);
                for(Node* n : poSeq){
                    int nIndex = n->getIndex();
                    if(n != activeT->getRoot()){
                        double v = activeT->getBranchLength(n);
                        transProb->setProbs(activeTP[nIndex], c, nIndex, activeSubspace, v);
                    }
                }
            });
        }

        executor.run(transProbTasks).wait();

        if(bufferSize < numCurrentCats + numAux){
            delete [] tempCLBuffer;
            delete [] tempRescaleBuffer;

            bufferSize = (2 * numAux) + numCurrentCats;
            tempCLBuffer = new double[numNodes * bufferSize * stateSpace];
            tempRescaleBuffer = new double[numNodes * bufferSize];
            std::fill(tempCLBuffer, tempCLBuffer + numNodes * bufferSize * stateSpace, 0.0);
            std::fill(tempRescaleBuffer, tempRescaleBuffer + numNodes * bufferSize, 0.0);
        }

        int nodeSpacer = bufferSize * stateSpace;
        int activeNodeSpacer = numTestableCats * stateSpace;
        int fullSpacer = numNodes * nodeSpacer;

        for(auto n : poSeq){
            int nIndex = n->getIndex();
            auto pNN = tempCLBuffer + nIndex * nodeSpacer;
            if(!n->getIsTip()){
                std::fill(pNN, pNN + activeNodeSpacer, 1.0);
                std::set<Node*>& nNeighbors = n->getNeighborRef();

                for(Node* d : nNeighbors){
                    if(d != n->getAncestor()){
                        int dIndex = d->getIndex();
                        double* pN = pNN;
                        double* pD = tempCLBuffer + (dIndex * nodeSpacer);

                        for(int c = 0; c < numTestableCats; c++){
                            Matrix<double> P = (*transProb)(activeTP[dIndex], c, dIndex);
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
                            pN+=(3-numClasses)*61;
                            pD+=stateSpace;
                        }
                    }
                }

                double* rescalePointer = tempRescaleBuffer + (nIndex * bufferSize);
                std::fill(rescalePointer, rescalePointer + bufferSize, 0.0);

                for(int c = 0; c < numTestableCats; c++){
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
            else{
                double* dataPointer = (*postOrder)(nIndex, activeCL[nIndex], 0) + randomSite*stateSpace;
                double* sitePointer = tempCLBuffer + nIndex*nodeSpacer;
                for(int c = 0; c < numTestableCats; c++){
                    for(int i = 0; i < stateSpace; i++){
                        *sitePointer = dataPointer[i];
                        sitePointer++;
                    }
                }
            }
        }


        int rIndex = activeT->getRoot()->getIndex();
        std::vector<double> f = rateMatrix->getStationary();
        std::vector<double> likelihoodVec(numTestableCats, 0.0);
        double* siteRoot = tempCLBuffer + (rIndex * nodeSpacer);

        for(int c = 0; c < numTestableCats; c++){
            for(int i = 0; i < activeSubspace; i++){
                likelihoodVec[c] += siteRoot[i]*f[i];
            }
            siteRoot += stateSpace;
            likelihoodVec[c] = std::log(likelihoodVec[c]);
        }

        double* rescaleBufferPointer = tempRescaleBuffer;
        for(int i = 0; i < numNodes; i++){
            for(int c = 0; c < numTestableCats; c++){
                likelihoodVec[c] += rescaleBufferPointer[c];
            }
            rescaleBufferPointer += bufferSize;
        }
        #else
        std::vector<double> likelihoodVec(numTestableCats, 0.0);
        #endif

        for(double& d : likelihoodVec){
            d *= samplingPower;
        }

        for(int c = 0; c < numCurrentCats; c++){
            likelihoodVec[c] += samplingPower * std::log(currentCategories[c].members.size());
        }

        for(int c = numCurrentCats; c < numTestableCats; c++){
            likelihoodVec[c] += samplingPower * alphaSplit;
        }

        double maxL = *std::max_element(likelihoodVec.begin(), likelihoodVec.end());
        double total = 0.0;
        for(double& d : likelihoodVec){
            d -= maxL;
            d = std::exp(d);
            total += d;
        }

        double categoryDraw = total * rng.uniformRv();
        bool assigned = false;

        total = 0.0;
        for(int i = 0; i < likelihoodVec.size(); i++){
            total += likelihoodVec[i];
            if(total > categoryDraw){
                if(i < numCurrentCats) { //It already exists
                    rateMatrix->assignMember(randomSite, i);
                }
                else { // Generate a new category
                    rateMatrix->addCategory(newOmegas[i - numCurrentCats]);
                    rateMatrix->assignMember(randomSite, numCurrentCats);
                    regenerateTransitionProbs(numCurrentCats);
                }
                assigned = true;
                break;
            }
        }

        if(assigned == false){
            std::cout << randomSite << " was not assigned to any category.\n";
            std::cout << "Log Likelihoods: \n";
            for(double& a : likelihoodVec){
                std::cout << a << "\n";
            }
            std::cout << std::flush;
            Msg::error("Failed to assign!");
        }

        int bufferDiff = transProb->getNumMatrices() - rateMatrix->getNumCategories();
        if(bufferDiff > 30){
            transProb->deleteNQ(bufferDiff - 30);
        }
    }

    rateMatrix->regenerateCatPrior();
    rateMatrix->regenerateAssignments();

    delete [] tempCLBuffer;
    delete [] tempRescaleBuffer;

    return INFINITY;
}

void DPCMMModel::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
}

void DPCMMModel::setCountTuningEvents(bool shouldCount) {
    tree->setCountTuningEvents(shouldCount);
    rateMatrix->setCountTuningEvents(shouldCount);
}


std::vector<double> DPCMMModel::getTunableParameterRecord() const {
    std::vector<double> record = {
        rateMatrix->getK(), rateMatrix->getR()
    };
    for(double entry : rateMatrix->getRawStationary())
        record.push_back(entry);
    for(double entry : tree->getTree()->getBranchProportions())
        record.push_back(entry);

    return record;
}

std::vector<double> DPCMMModel::getTunableParameters() const {
    std::vector<double> returnVec(5, 0.0);
    returnVec[0] = tree->getBranchAlpha();
    returnVec[1] = tree->getTreeDelta();
    returnVec[2] = rateMatrix->getStationaryAlpha();
    returnVec[3] = rateMatrix->getKDelta();
    returnVec[4] = rateMatrix->getRDelta();
    return returnVec;
}

void DPCMMModel::setTunableParameters(const std::vector<double> & v){
    tree->setBranchAlpha(v[0]);
    tree->setTreeDelta(v[1]);
    rateMatrix->setStationaryAlpha(v[2]);
    rateMatrix->setKDelta(v[3]);
    rateMatrix->setRDelta(v[4]);
}

void DPCMMModel::printAcceptanceRates(){
    std::cout << "Tree Acceptance Rate: " << tree->getTreeRate() << "\tBranch Acceptance Rate: " << tree->getBranchRate() <<
    "\tStationary Acceptance Rate: " << rateMatrix->getStationaryRate() << "\tK Acceptance Rate: " << rateMatrix->getKRate() <<
    "\tOmega Acceptance Rate: " << rateMatrix->getOmegaRate() << "\tR Acceptance Rate: " << rateMatrix->getRRate() << std::endl;
}

void DPCMMModel::printTabular(int i){
    if(i == 0){
        std::string returnString = "Iteration\tPosterior\tLikelihood\tPrior\tTreeLength\tOmegaCount\tK\tR";
        for(int i = 0; i < 61; i++)
            returnString += "\tPi[" + std::to_string(i) + "]";
            

        std::cout << returnString << std::endl;
    }
    else{
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                                std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" + std::to_string(tree->getTree()->getTreeLength()) + "\t" +
                                std::to_string(rateMatrix->getActiveOmegas()) + "\t" +
                                std::to_string(rateMatrix->getK()) + "\t" +std::to_string(rateMatrix->getR());
        std::vector<double> stationary = rateMatrix->getRawStationary();
        for(double i : stationary){
            returnString += "\t" + std::to_string(i);
        }

        std::cout << returnString << std::endl;
    }
}

void DPCMMModel::writeLogHeaders(){
    if(analysisLog != ""){
        std::string tabHeader = "Iteration\tPosterior\tLikelihood\tPrior\tTreeLength\tOmegaCount\tK\tR";
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

    if(dppLog != ""){
        std::string dppHeader = "Iteration\tPosterior\tCategoryCount";
        for(int i = 0; i < numChar; i++)
            dppHeader += "\tOmega[" + std::to_string(i) + "]" + "\tOmegaIncrement1[" + std::to_string(i) + "]\tOmegaIncrement2[" + std::to_string(i) + "]";
        dppHeader += "\n";

        std::ofstream outFile(dppLog, std::ios::out);
        outFile << dppHeader;
    }

    #ifndef SAMPLE_PRIOR
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
    #endif
}

void DPCMMModel::writeLogData(int i) {
    if(analysisLog != ""){
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                                std::to_string(currentLikelihood) + "\t" + std::to_string(lnPrior()) + "\t" + std::to_string(tree->getTree()->getTreeLength()) + "\t" +
                                std::to_string(rateMatrix->getActiveOmegas()) + "\t" +
                                std::to_string(rateMatrix->getK()) + "\t" + std::to_string(rateMatrix->getR());
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

    if(dppLog != ""){
        std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t";
        std::vector<Category> categories = rateMatrix->getCategories();
        returnString += std::to_string(categories.size());
        std::vector<int> assignments = rateMatrix->getAssignments();
        for(int c : assignments){
            returnString += "\t" + std::to_string(categories[c].omegas[0]) + "\t" + std::to_string(categories[c].omegas[1]) + "\t" + std::to_string(categories[c].omegas[2]);
        }
        returnString += "\n";

        std::ofstream outFile(dppLog, std::ios::app);
        outFile << returnString;
    }

    #ifndef SAMPLE_PRIOR
    if(tipsLog != "" || ancestralLog != ""){
        RandomVariable& rng = RandomVariable::randomVariableInstance();

        std::string tipString = std::to_string(i);
        std::string ancestralString = std::to_string(i);
        std::vector<Node*> tips = tree->getTree()->getTips();

        TreeObject* activeT = tree->getTree();
        std::vector<Node*> preOrderSeq = activeT->getPostOrderSeq();
        std::reverse(preOrderSeq.begin(), preOrderSeq.end());
        std::vector<int> assignments = rateMatrix->getAssignments();
        std::vector<Category> categories = rateMatrix->getCategories();
        int numCats = rateMatrix->getNumCategories();

        std::vector<std::array<double, 3>> dNdSVec;
        for(int i = 0; i < numCats; i++){
            std::array<double, 3> dNdS = rateMatrix->dNdS(i);
            dNdSVec.push_back(dNdS);
        }

        int numClasses = rateMatrix->getActiveOmegas();

        int* reconstructedStates = new int[numNodes*numChar];
        double* reconstructeddNdS = new double[numNodes*numChar];

        std::fill(reconstructeddNdS, reconstructeddNdS + numNodes*numChar, 0.0);
        std::fill(reconstructedStates, reconstructedStates + numNodes*numChar, -1);

        tf::Taskflow phyloTaskflow;
        std::unordered_map<int, tf::Task> taskMap;

        Node* root = activeT->getRoot();
        taskMap.insert(std::make_pair(root->getIndex(), phyloTaskflow.emplace([this, root, &rng, reconstructedStates, reconstructeddNdS, assignments, dNdSVec, numClasses](){
            int rIndex = root->getIndex();
            double* pR = (*postOrder)(rIndex, activeCL[rIndex], 0);
            int* reconstructedP = reconstructedStates + rIndex*numChar;
            double* dNdSP = reconstructeddNdS + rIndex*numChar;

            for(int c = 0; c < numChar; c++){
                double total = 0;
                for(int i = 0; i < numClasses * 61; i++){
                    total += pR[i];
                }

                double draw = rng.uniformRv() * total;

                double sum = 0;
                bool success = false;
                for(int i = 0; i < numClasses * 61; i++){
                    sum += pR[i];

                    if(sum >= draw){
                        *reconstructedP = i;
                        *dNdSP = dNdSVec[assignments[c]][(int)(i/61.0)];
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
                taskMap.insert(std::make_pair(n->getIndex(), phyloTaskflow.emplace([this, n, &rng, reconstructedStates, reconstructeddNdS, assignments, dNdSVec, numClasses](){
                    int nIndex = n->getIndex();
                    double* pN = (*postOrder)(nIndex, activeCL[nIndex], 0);

                    int* reconstructedP = reconstructedStates + nIndex*numChar;
                    double* dNdSP = reconstructeddNdS + nIndex*numChar;

                    int ancestorIndex = n->getAncestor()->getIndex();

                    for(int c = 0; c < numChar; c++){
                        Matrix<double> P = (*transProb)(activeTP[nIndex], assignments[c], nIndex);
                        int ancestorState = *(reconstructedStates + ancestorIndex*numChar + c);

                        double total = 0;
                        for(int i = 0; i < numClasses*61; i++){
                            total += P(ancestorState, i) * pN[i];
                        }

                        double draw = rng.uniformRv() * total;

                        double sum = 0;
                        bool success = false;
                        for(int i = 0; i < numClasses*61; i++){
                            sum += P(ancestorState, i) * pN[i];

                            if(sum >= draw){
                                *reconstructedP = i;
                                *dNdSP = dNdSVec[assignments[c]][(int)(i/61.0)];
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
            outFile << tipString << std::endl;
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
            outFile << ancestralString << std::endl;
        }

        delete [] reconstructedStates;
        delete [] reconstructeddNdS;
    }
    #endif
}
