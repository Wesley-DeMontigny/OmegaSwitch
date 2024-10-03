#include "ConditionalLikelihood.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"

ConditionalLikelihood::ConditionalLikelihood(Alignment* aln, int nR) : numNodes(aln->getNumTaxa() * 2), numRates(nR) {
    numChar = aln->getNumChar();
    stateSpace = aln->getStateSpace();
    int width = numNodes*numChar*stateSpace*numRates;
    condLikelihoods[0] = new double[2 * width];
    condLikelihoods[1] = condLikelihoods[0] + (width);

    activeCLs = new int[numNodes];
    for(int i = 0; i < numNodes; i++)
        activeCLs[i] = 0;

    for(int i = 0; i < width; i++){
        condLikelihoods[0][i] = 0.0;
        condLikelihoods[1][i] = 0.0;
    }

    for(int index = 0; index < aln->getNumTaxa(); index++){
        for(int r = 0; r < numRates; r++){
            double* p = (*this)(index, r, 0);
            for(int i = 0; i < numChar; i++){
                unsigned long long int state = aln->getMatrix()[index][i];

                unsigned long long int mask = 1;
                bool assigned = false;
                for(int j = 0; j < stateSpace; j++) {
                    if((mask & state) != 0){
                        *p = 1.0;
                        assigned = true;
                    }
                    mask <<= 1;
                    p++;
                }

                if(assigned == false){
                    Msg::error("Never assigned a conditional value at (" + std::to_string(index) + ", " + std::to_string(i) + ")! This has state value " + std::to_string(state));
                }
            }
        }
    }
}

ConditionalLikelihood::ConditionalLikelihood(int nT, int nC, int nR, int s) : numNodes(nT * 2), stateSpace(s), numRates(nR), numChar(nC) {
    int width = numNodes*numChar*stateSpace*numRates;
    condLikelihoods[0] = new double[2 * width];
    condLikelihoods[1] = condLikelihoods[0] + (width);

    activeCLs = new int[numNodes];
    for(int i = 0; i < numNodes; i++)
        activeCLs[i] = 0;

    for(int i = 0; i < width; i++){
        condLikelihoods[0][i] = 0.0;
        condLikelihoods[1][i] = 0.0;
    }
}

ConditionalLikelihood::~ConditionalLikelihood(){
    delete [] condLikelihoods[0];
    delete activeCLs;
}


double* ConditionalLikelihood::operator()(int n, int s, int r){
    return condLikelihoods[s] + n*numChar*stateSpace + (r*numNodes*numChar*stateSpace);
}

double* ConditionalLikelihood::operator[](int n){
    return condLikelihoods[activeCLs[n]] + n*numChar*stateSpace;
}

void ConditionalLikelihood::flipCL(int n){
    activeCLs[n] ^= 1;
}