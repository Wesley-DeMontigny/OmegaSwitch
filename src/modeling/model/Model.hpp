#ifndef PHYLO_CTMC_HPP
#define PHYLO_CTMC_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class Alignment;
class RandomVariable;
class RateMatrix;
class CodonMultiMatrix;
class DirichletProcessPrior;

class Model {
    public:
        Model(void) = delete;
        Model(Alignment* a, TreeParameter* t, CodonMultiMatrix* m, DirichletProcessPrior* d);
        ~Model();
        double lnLikelihood() {return currentLikelihood;}
        double lnPrior();
        double regenerateIntoLikelihoodBuffer(int site, int category, bool update);
        void forceRegenerate(int site, int category, bool update);
        void regenerateLikelihood();
        TransitionProbability* getTransitionProbability() { return transProb; }
        void accept();
        void reject();
        void tuneMoves();
        std::string tabularOut(int i);
        std::string tabularHeader();
        std::string treeOut(int i);
        std::string treeHeader();
        std::string dppOut(int i);
        std::string dppHeader();
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
        void postOrderPrune();
};

#endif