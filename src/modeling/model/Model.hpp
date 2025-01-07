#ifndef PHYLO_CTMC_HPP
#define PHYLO_CTMC_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "core/Alignment.hpp"

class ConditionalLikelihood;
class TransitionProbability;
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

        void regenerateLikelihood();
        void regenerateTransitionProbs(int site, int category);
        void reconstructTips();

        double testCategory(int site, int category, bool update);

        int getNumTaxa(){return aln->getNumTaxa();}
        int getNumChar(){return numChar;}
        int getNumNodes(){return numNodes;}

        TransitionProbability* getTransitionProbability() { return transProb; }
        ConditionalLikelihood* getConditionalLikelihood() { return postOrder; }

        void accept();
        void reject();
        void tuneMoves();

        std::string tabularHeader();
        std::string tabularOut(int i);
        std::string treeHeader();
        std::string treeOut(int i);
        std::string dppHeader();
        std::string dppOut(int i);
        std::string tipsHeader();
        std::string tipsOut(int i);
    protected:
        double oldLikelihood;
        double currentLikelihood;
    private:
        int stateSpace;
        int numChar;
        int numNodes;
        bool* activeTP;
        bool* activeCL;
        tf::Executor executor;
        CodonMultiMatrix* rateMatrix;
        Alignment* aln;
        ConditionalLikelihood* postOrder;
        TransitionProbability* transProb;
        double* reconstruction;
        double* rescaling;
        TreeParameter* tree;
        DirichletProcessPrior* dpp;
};

#endif