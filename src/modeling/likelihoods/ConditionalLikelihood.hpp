#ifndef CONDITIONAL_LIKELIHOOD_HPP
#define CONDITIONAL_LIKELIHOOD_HPP
#include <iostream>

class Alignment;
class Node;

class ConditionalLikelihood{
    public:
        ConditionalLikelihood(void) = delete;
        ConditionalLikelihood(Alignment* aln, int nR);
        ConditionalLikelihood(int nT, int nC, int nR, int s);
        ~ConditionalLikelihood();
        double* operator()(int n, int s, int r);
        double* operator[](int n);
        void flipCL(int n);
    private:
        double* condLikelihoods[2];
        int* activeCLs;
        int numChar;
        int numNodes;
        int numRates;
        int stateSpace;
};

#endif