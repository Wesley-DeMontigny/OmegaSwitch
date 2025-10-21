#ifndef RJ_DPP_MODEL_HPP
#define RJ_DPP_MODEL_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "core/Alignment.hpp"
#include "Model.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class RJDPPMatrix;
class Settings;

class RJDPPModel : public Model {
    public:
        RJDPPModel(void) = delete;
        RJDPPModel(Settings* s, Alignment* a, TreeParameter* t, RJDPPMatrix* m, tf::Executor& e);
        ~RJDPPModel();
        
        double                  lnLikelihood() override {return currentLikelihood; }                                    //
        double                  lnPrior() override;                                                                     //
        double updateDPP();
        std::vector<double>     getTunableParameterRecord() override;                                                   //
        std::vector<double>     getTunableParameters() override;                                                        //
        void                    accept() override;                                                                      //
        void                    printAcceptanceRates() override;                                                        //
        void                    printTabular(int i) override;                                                           //
        void                    regenerateLikelihood() override;                                                        //
        void                    reject() override;                                                                      //
        void                    setTunableParameters(const std::vector<double> & v) override;                           //
        void                    tuneMoves() override;                                                                   //
        void                    writeLogData(int i) override;                                                           //
        void                    writeLogHeaders() override;                                                             //
        void                    regenerateTransitionProbs(int category);
    private:
        Alignment* aln;
        bool* activeCL;
        bool* activeTP;
        ConditionalLikelihood* postOrder;
        double currentLikelihood;
        double oldLikelihood;
        double omegaLambda;
        double* rescaling;
        int numChar;
        int numGibbsUpdate;
        int numNodes;
        int stateSpace;
        RJDPPMatrix* rateMatrix;
        std::string             analysisLog;
        std::string             ancestralLog;
        std::string             branchLog;
        std::string             dppLog;
        std::string             tipsLog;
        std::string             treeLog;
        tf::Executor& executor;
        TransitionProbability* transProb;
        TreeParameter* tree;
};

#endif