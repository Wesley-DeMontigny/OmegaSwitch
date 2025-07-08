#ifndef RJ_MODEL_HPP
#define RJ_MODEL_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "core/Alignment.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class RJMatrix;
class Settings;

class RJModel {
    public:
        RJModel(void) = delete;
        RJModel(Settings s, Alignment* a, TreeParameter* t, RJMatrix* m);
        ~RJModel();

        void accept();
        void reject();
        void tuneMoves();
        double lnLikelihood() {return currentLikelihood;}
        double lnPrior();

        void regenerateLikelihood();

        int getNumTaxa(){return aln->getNumTaxa();}
        int getNumChar(){return numChar;}
        int getNumNodes(){return numNodes;}

        TransitionProbability* getTransitionProbability() { return transProb; }
        ConditionalLikelihood* getConditionalLikelihood() { return postOrder; }

        std::string tabularHeader();
        std::string tabularOut(int i);
        std::string treeHeader();
        std::string treeOut(int i);
        std::string tipsHeader();
        std::string ancestralHeader();
        std::tuple<std::string, std::string> reconstructionOut(int i);
        std::string branchHeader();
        std::string branchOut(int i);
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

        RJMatrix* rateMatrix;
        Alignment* aln;
        ConditionalLikelihood* postOrder;
        TransitionProbability* transProb;
        TreeParameter* tree;
};

#endif