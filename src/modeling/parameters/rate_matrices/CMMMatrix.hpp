#ifndef CMM_MATRIX_HPP
#define CMM_MATRIX_HPP
#include "misc/Matrix.hpp"
#include "misc/Msg.hpp"
#include "MatrixHelper.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <memory>
#include <array>
#include <set>
#include <vector>

class Settings;

using namespace MatrixHelper;

/**
 * @brief 
 */
class CMMMatrix : public Parameter {
    public:
                                    CMMMatrix(Settings& settings);                           //

        double                      getActiveOmegas() const {return currentActiveOmegas;}   //
        double                      getK() const {return currentParams[0];}                 //
        double                      getKDelta() const {return tuningState->kDelta;}         //
        double                      getKRate() const;                                       //
        double                      getOmega(int i) const {return currentParams[2+i];}      // Get the ith Omega (or Omega increment)
        double                      getOmegaDelta() const {return tuningState->omegaDelta;} //
        double                      getOmegaRate() const;                                   //
        double                      getR() const {return currentParams[1];}                 //
        double                      getRDelta() const {return tuningState->rDelta;}         //
        double                      getRRate() const;                                       //
        bool                        isRegimeCountFixed() const {return fixedRegimes > 0;}   //
        double                      getStationaryAlpha() const {return tuningState->stationaryAlpha;}//
        double                      getStationaryRate() const;                              //
        double                      lnPrior() override;                                     //
        double                      updateActiveOmegas();                                   //
        double                      updateK();                                              //
        double                      updateOmega();                                          //
        double                      updateR();                                              //
        double                      updateStationary();                                     //
        Matrix<double>              Q() const;                                              //
        std::array<double, 5>       dNdS() const;                                           //
        const std::vector<double>&  getRawStationary() const {return currentStationary;}    //
        std::vector<double>         getStationary() const;                                  //
        void                        accept() override;                                      //
        void                        reject() override;                                      //
        void                        setCountTuningEvents(bool shouldCount) {countTuningEvents = shouldCount;}//
        void                        setActiveOmegas(int o) {                                //
                                    currentActiveOmegas = o; 
                                    rebuildQMatrix();
                                    }
        void                        setKDelta(double delta) {tuningState->kDelta = delta;}  //
        void                        setOmegaDelta(double delta) {tuningState->omegaDelta = delta;}//
        void                        setRDelta(double delta) {tuningState->rDelta = delta;}  //
        void                        setStationaryAlpha(double alpha) {tuningState->stationaryAlpha = alpha;}//
        void                        shareTuningWith(CMMMatrix& m) {tuningState = m.tuningState;}//
        void                        tune() override;                                        //
    private:

        struct CMMTuningState {
            double kDelta = 0.5;
            double omegaDelta = 0.5;
            double rDelta = 0.5;
            double stationaryAlpha = 1000.0;
            ProposalTuningStats kStats;
            ProposalTuningStats omegaStats;
            ProposalTuningStats rStats;
            ProposalTuningStats stationaryStats;
        };

        void                        rebuildQMatrix();                                       //
        double                      currentStationaryPrior = 0;                             //
        const double                kLambda;                                                //
        double                      oldStationaryPrior = 0;                                 //
        const double                omegaLambda;                                            //
        const double                rLambda;                                                //
        const int                   fixedRegimes;                                           //
        int                         currentActiveOmegas = 1;                                //
        int                         kAcceptCount = 0;                                       //
        int                         kCount = 0;                                             //
        bool                        countTuningEvents = true;                               //
        int                         oldActiveOmegas = 1;                                    //
        int                         omegaAcceptCount = 0;                                   //
        int                         omegaCount = 0;                                         //
        int                         rAcceptCount = 0;                                       //
        int                         rCount = 0;                                             //
        int                         stationaryAcceptCount = 0;                              //
        int                         stationaryCount = 0;                                    //
        Matrix<double>              currentQMatrix;                                         //
        Matrix<double>              oldQMatrix;                                             //
        MatrixMoves                 moveChoice = MatrixMoves::NO_MOVE;                      //
        std::array<double, 7>       currentParamPriors = {0,0,0,0,0,0,0};                   //
        std::array<double, 7>       currentParams = {0,0,0,0,0,0,0};                        // Parameters for the rate matrix in the order of K, R, O1, O2, O3, O4, O5
        std::array<double, 7>       oldParamPriors = {0,0,0,0,0,0,0};                       // Old rate matrix parameters
        std::array<double, 7>       oldParams = {0,0,0,0,0,0,0};                            //
        std::vector<double>         currentStationary;                                      //
        std::vector<double>         oldStationary;                                          //
        std::vector<double>         stationaryPriorAlpha;                                   //
        std::vector<int>            randomStates;                                           //
        std::shared_ptr<CMMTuningState> tuningState = std::make_shared<CMMTuningState>();  //
};

#endif
