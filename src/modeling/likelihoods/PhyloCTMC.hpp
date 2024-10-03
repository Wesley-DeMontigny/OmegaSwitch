#ifndef PHYLO_CTMC_HPP
#define PHYLO_CTMC_HPP
#include "modeling/likelihoods/LikelihoodNode.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class Alignment;
class RandomVariable;
class RateMatrix;
class CodonMultiMatrix;
class DirichletProcessPrior;

class PhyloCTMC : public LikelihoodNode{
    public:
        PhyloCTMC(void) = delete;
        PhyloCTMC(Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d);
        ~PhyloCTMC();
        double lnLikelihood() {return currentLikelihood;}
        double regenerateAtSite(int site, int category, bool update);
        void regenerate();
        void accept();
        void reject();
        TransitionProbability* getTransitionProbability() {return transProb;}
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
        TransitionProbability* transProb;
        TreeParameter* tree;
        DirichletProcessPrior* dpp;
};

#endif