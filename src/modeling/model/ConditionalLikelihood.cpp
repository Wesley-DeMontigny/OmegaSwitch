#include "ConditionalLikelihood.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"

ConditionalLikelihood::ConditionalLikelihood(Alignment* aln, int nN, int nR, int ss) : numNodes(nN), numRates(nR), stateSpace(ss) {
    numChar = aln->getNumChar();
    int width = numNodes*numChar*stateSpace*numRates;
    condLikelihoods[0] = new double[2 * width];
    condLikelihoods[1] = condLikelihoods[0] + (width);

    for(int i = 0; i < width; i++){
        condLikelihoods[0][i] = 0.0;
        condLikelihoods[1][i] = 0.0;
    }

    int unseenRates = ss / 61;

    for(int index = 0; index < aln->getNumTaxa(); index++){
        for(int r = 0; r < numRates; r++){
            double* p = (*this)(index, 0, r);
            for(int i = 0; i < numChar; i++){
                const std::bitset<61>& state = aln->getMatrix()[index][i];

                bool assigned = false;
                for(int j = 0; j < 61; j++) {
                    if(state[j] == 1){
                        for(int u = 0; u < unseenRates; u++){
                            *(p + (61 * u)) = 1.0;
                        }
                        assigned = true;
                    }
                    p++;
                }
                p += 61 * (unseenRates-1);

                if(assigned == false){
                    Msg::error("Never assigned a conditional value at (" + std::to_string(index) + ", " + std::to_string(i) + ")!");
                }
            }
        }
    }
}

ConditionalLikelihood::~ConditionalLikelihood(){
    delete [] condLikelihoods[0];
}


double* ConditionalLikelihood::operator()(int n, int s, int r){
    return condLikelihoods[s] + n*numChar*stateSpace + (r*numNodes*numChar*stateSpace);
}
