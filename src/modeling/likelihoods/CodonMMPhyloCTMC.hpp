#ifndef CODOM_MM_PHYLO_CTMC_HPP
#define CODOM_MM_PHYLO_CTMC_HPP
#include "modeling/likelihoods/LikelihoodNode.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"

class ConditionalLikelihood;
class MultiMatrixTransitionProbability;
class Alignment;
class RandomVariable;
class RateMatrix;
class CodonMultiMatrix;
class DirichletProcessPrior;

class CodonMMPhyloCTMC : public LikelihoodNode{
    public:
        CodonMMPhyloCTMC(void) = delete;
        CodonMMPhyloCTMC(Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d);
        ~CodonMMPhyloCTMC();
        double lnLikelihood() {return currentLikelihood;}
        double regenerateAtSite(int site, int category, bool update);
        void regenerate();
        void accept();
        void reject();
        MultiMatrixTransitionProbability* getTransitionProbability() {return transProb;}
        ConditionalLikelihood* getPostOrderL() {return postOrder;}
        std::string writeValue() {return std::to_string(currentLikelihood);}
    protected:
        double oldLikelihood;
        double currentLikelihood;
    private:
        int stateSpace;
        bool* activeTP;
        bool* activeCL;
        CodonMultiMatrix* rateMatrix;
        Alignment* aln;
        ConditionalLikelihood* postOrder;
        MultiMatrixTransitionProbability* transProb;
        TreeParameter* tree;
        DirichletProcessPrior* dpp;
};

#endif