#ifndef M1_Model_HPP
#define M1_Model_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "core/Alignment.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class M1Matrix;
class Settings;

class M1Model {
    public:
        M1Model(void) = delete;
        M1Model(Settings s, Alignment* a, TreeParameter* t, M1Matrix* m);
        ~M1Model();

        double lnLikelihood() {return currentLikelihood;}
        double lnPrior();

        void regenerateLikelihood();

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
    protected:
        double oldLikelihood;
        double currentLikelihood;
    private:
        int stateSpace;
        int numChar;
        int numNodes;
        bool* activeTP;
        bool* activeCL;
        double* rescaling;

        tf::Executor executor;

        M1Matrix* rateMatrix;
        Alignment* aln;
        ConditionalLikelihood* postOrder;
        TransitionProbability* transProb;
        TreeParameter* tree;
};

#endif