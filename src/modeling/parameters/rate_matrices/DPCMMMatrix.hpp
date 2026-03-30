#ifndef DP_CMM_MATRIX_HPP
#define DP_CMM_MATRIX_HPP
#include "misc/Matrix.hpp"
#include "misc/Msg.hpp"
#include "misc/Probability.hpp"
#include "MatrixHelper.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <memory>
#include <array>
#include <optional>
#include <set>
#include <vector>

class Settings;

struct Category {
    std::array<double, 3> omegas;
    std::vector<int> members;
    bool dirty;
};

using namespace MatrixHelper;

/**
 * @brief 
 * 
 */
class DPCMMMatrix : public Parameter {
    public:
                                        DPCMMMatrix(Settings& settings, int nC);                        //

        double                          getK() const {return currentParams[0];}                         //
        double                          getKDelta() const {return tuningState->kDelta;}                 //
        double                          getKRate() const;                                               //
        double                          getOmegaDelta() const {return tuningState->omegaDelta;}         //
        double                          getOmegaRate() const;                                           //
        double                          getR() const {return currentParams[1];}                         //
        double                          getRDelta() const {return tuningState->rDelta;}                 //
        double                          getRRate()const ;                                               //
        bool                            isRegimeCountFixed() const {return fixedRegimes > 0;}          //
        double                          getStationaryAlpha() const {return tuningState->stationaryAlpha;}// 
        double                          getStationaryRate() const;                                      //
        double                          lnPrior() override;                                             //
        double                          updateActiveOmegas();                                           //
        double                          updateK();                                                      //
        double                          updateOmega();                                                  //
        double                          updateR();                                                      //
        double                          updateStationary();                                             //
        double                          getAlpha() const {return dpAlpha;}                              //
        double                          getActiveOmegas() const {return currentActiveOmegas;}           //
        std::optional<int>              unassignMember(int member);                                     //
        int                             getNumCategories() const {return currentCategories.size();}     //
        Matrix<double>                  Q(int i) const;                                                 //
        Matrix<double>                  Q(const std::array<double, 3>& omegas) const;                   //
        std::array<double, 3>           dNdS(int i) const;                                              //
        const std::vector<Category>&    getCategories() const {return currentCategories;}               //
        const std::vector<double>&      getRawStationary() const {return currentStationary;}            //
        std::vector<double>             getStationary(int omegaCount);                                  //
        const std::vector<int>&         getAssignments() const { return assignments;}                   //
        void                            accept() override;                                              //
        void                            reject() override;                                              //
        void                            setCountTuningEvents(bool shouldCount) {countTuningEvents = shouldCount;}//
        void                            setActiveOmegas(int o) {                                        //
                                        currentActiveOmegas = o; 
                                        rebuildQMatrix();
                                        }
        void                            setKDelta(double delta) {tuningState->kDelta = delta;}         //
        void                            setOmegaDelta(double delta) {tuningState->omegaDelta = delta;} //
        void                            setRDelta(double delta) {tuningState->rDelta = delta;}         //
        void                            setStationaryAlpha(double alpha) {tuningState->stationaryAlpha = alpha;}//
        void                            shareTuningWith(DPCMMMatrix& m) {tuningState = m.tuningState;} //
        void                            tune() override;                                                //
        void                            addCategory(const std::array<double, 3>& omegas);               //
        void                            assignMember(int member, int category);                         //
        void                            popCategories(int n);                                           //
        void                            regenerateAssignments();                                        //
        void                            removeCategory(int index);                                      //
        void                            regenerateCatPrior();                                           //
    private:
        struct ProposalTuningStats {
            int acceptCount = 0;
            int count = 0;
        };
        
        struct DPCMMTuningState {
            double kDelta = 0.5;
            double rDelta = 0.5;
            double omegaDelta = 0.5;
            double stationaryAlpha = 1000.0;
            ProposalTuningStats kStats;
            ProposalTuningStats rStats;
            ProposalTuningStats omegaStats;
            ProposalTuningStats stationaryStats;
        };

        void                            rebuildQMatrix();                                               //
        double                          calculateAlpha(double expectedCat, int members) const;          //
        double                          expectedCategories(double a, int members) const;                //

        double                          currentCatPrior = 0.0;                                          //
        double                          currentStationaryPrior = 0;                                     //
        double                          dpAlpha = 0.0;                                                  //
        const double                    kLambda;                                                        //
        double                          oldCatPrior = 0.0;                                              //
        double                          oldStationaryPrior = 0;                                         //
        const double                    omegaLambda;                                                    //
        const double                    rLambda;                                                        //
        const int                       fixedRegimes;                                                   //
        id_t                            omegaCount = 0;                                                 //
        int                             currentActiveOmegas = 1;                                        //
        int                             kAcceptCount = 0;                                               //
        int                             kCount = 0;                                                     //
        bool                            countTuningEvents = true;                                       //
        int                             numChar;                                                        //
        int                             oldActiveOmegas = 1;                                            //
        int                             omegaAcceptCount = 0;                                           //
        int                             rAcceptCount = 0;                                               //
        int                             rCount = 0;                                                     //
        int                             stationaryAcceptCount = 0;                                      //
        int                             stationaryCount = 0;                                            //
        Matrix<double>                  currentQMatrix;                                                 //
        Matrix<double>                  oldQMatrix;                                                     //
        MatrixMoves                     moveChoice = MatrixMoves::NO_MOVE;                              //
        std::array<double, 2>           currentParamPriors = {0,0};                                     //
        std::array<double, 2>           currentParams = {0,0};                                          // Parameters for the rate matrix in the order of K, R
        std::array<double, 2>           oldParamPriors = {0,0};                                         //
        std::array<double, 2>           oldParams = {0,0};                                              //
        std::vector<Category>           currentCategories;                                              //
        std::vector<Category>           oldCategories;                                                  //
        std::vector<double>             currentStationary;                                              //
        std::vector<double>             oldStationary;                                                  //
        std::vector<double>             stationaryPriorAlpha;                                           //
        std::vector<int>                assignments;                                                    //
        std::vector<int>                randomStates;                                                   //
        std::shared_ptr<DPCMMTuningState> tuningState = std::make_shared<DPCMMTuningState>();          //
};

#endif
