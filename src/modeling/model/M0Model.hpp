#ifndef M0_MODEL_HPP
#define M0_MODEL_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "misc/Alignment.hpp"
#include "Model.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class M0Matrix;
class Settings;

/**
 * @brief 
 * 
 */
class M0Model : public Model {
    public:
                                M0Model(void) = delete;
                                M0Model(Settings* s, Alignment* a, TreeParameter* t, M0Matrix* m, tf::Executor& e);
                                ~M0Model();

        double                  lnLikelihood() override {return currentLikelihood; }                                    //
        double                  lnPrior() override;                                                                     //
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
    private:
        Alignment*              aln;
        bool*                   activeCL;
        bool*                   activeTP;
        ConditionalLikelihood*  postOrder;
        double                  currentLikelihood;
        double                  oldLikelihood;
        double*                 rescaling;
        int                     numChar;
        int                     numNodes;
        int                     stateSpace;
        M0Matrix*               rateMatrix;
        std::string             analysisLog;
        std::string             treeLog;
        tf::Executor&           executor;
        TransitionProbability*  transProb;
        TreeParameter*          tree;
};

#endif
