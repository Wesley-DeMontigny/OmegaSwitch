#ifndef RJ_MODEL_HPP
#define RJ_MODEL_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "core/Alignment.hpp"
#include "Model.hpp"
#include <taskflow/taskflow.hpp>
#include <string>

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class RJMatrix;
class Settings;

/**
 * @brief 
 * 
 */
class RJModel : public Model {
    public:
                                RJModel(void) = delete;
                                RJModel(Settings* s, Alignment* a, TreeParameter* t, RJMatrix* m, tf::Executor& e);     //
                                ~RJModel();                                                                             //

        double                  lnLikelihood() override {return currentLikelihood; }                                    //
        double                  lnPrior() override;                                                                     //
        std::vector<double>     getTunableParameterRecord() const override;                                             //
        std::vector<double>     getTunableParameters() const override;                                                  //
        void                    accept() override;                                                                      //
        void                    printAcceptanceRates() override;                                                        //
        void                    printTabular(int i) override;                                                           //
        void                    regenerateLikelihood() override;                                                        //
        void                    reject() override;                                                                      //
        void                    setTunableParameters(const std::vector<double> & v) override;                           //
        void                    tuneMoves() override;                                                                   //
        void                    writeLogData(int i) override;                                                           //
        void                    writeLogHeaders() override;                                                             //
    private:
        Alignment*              aln;                                                                                    //
        bool*                   activeCL;                                                                               //
        bool*                   activeTP;                                                                               //
        ConditionalLikelihood*  postOrder;                                                                              //
        double                  currentLikelihood;                                                                      //
        double                  oldLikelihood;                                                                          //
        double*                 rescaling;                                                                              //
        int                     numChar;                                                                                //
        int                     numNodes;                                                                               //
        int                     stateSpace;                                                                             //
        RJMatrix*               rateMatrix;                                                                             //
        std::string             analysisLog;                                                                            //
        std::string             treeLog;                                                                                //
        std::string             tipsLog;                                                                                //
        std::string             ancestralLog;                                                                           //
        std::string             branchLog;                                                                              //
        tf::Executor&           executor;                                                                               //
        TransitionProbability*  transProb;                                                                              //
        TreeParameter*          tree;                                                                                   //
};

#endif