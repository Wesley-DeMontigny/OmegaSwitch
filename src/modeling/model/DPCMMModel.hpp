#ifndef DP_CMM_MODEL_HPP
#define DP_CMM_MODEL_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "misc/Alignment.hpp"
#include "Model.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class DPCMMMatrix;
class Settings;

/**
 * @brief 
 * 
 */
class DPCMMModel : public Model {
    public:
        DPCMMModel(void) = delete;
        DPCMMModel(Settings* s, Alignment* a, TreeParameter* t, DPCMMMatrix* m, tf::Executor& e);
        ~DPCMMModel();
        
        double                  lnLikelihood() override {return currentLikelihood; }                                    //
        double                  lnPrior() override;                                                                     //
        double updateDPP();
        std::vector<double>     getTunableParameterRecord() const override;                                             //
        std::vector<double>     getTunableParameters() const override;                                                  //
        void                    accept() override;                                                                      //
        void                    printAcceptanceRates() override;                                                        //
        void                    printTabular(int i) override;                                                           //
        void                    regenerateLikelihood() override;                                                        //
        void                    reject() override;                                                                      //
        void                    setCountTuningEvents(bool shouldCount) override;                                        //
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
        DPCMMMatrix* rateMatrix;
        std::string             analysisLog;
        std::string             ancestralLog;
        std::string             dppLog;
        std::string             tipsLog;
        std::string             treeLog;
        tf::Executor& executor;
        TransitionProbability* transProb;
        TreeParameter* tree;
};

#endif
