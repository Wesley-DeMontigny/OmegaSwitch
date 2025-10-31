#ifndef RJ_MATRIX_HPP
#define RJ_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "MatrixHelper.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <array>
#include <set>
#include <vector>

class Settings;

using namespace MatrixHelper;

/**
 * @brief 
 */
class RJMatrix : public Parameter {
    public:
                                    RJMatrix(Settings& settings);                           //

        double                      getActiveOmegas() const {return currentActiveOmegas;}   //
        double                      getK() const {return currentParams[0];}                 //
        double                      getKRate() const;                                       //
        double                      getOmega(int i) const {return currentParams[2+i];}      // Get the ith Omega (or Omega increment)
        double                      getOmegaRate() const;                                   //
        double                      getR() const {return currentParams[1];}                 //
        double                      getRRate() const;                                       //
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
        void                        setActiveOmegas(int o) {                                //
                                    currentActiveOmegas = o; 
                                    rebuildQMatrix();
                                    }
        void                        tune() override;                                        //

        double                      kDelta;                                                 //
        double                      omegaDelta;                                             //
        double                      rDelta;                                                 //
        double                      stationaryAlpha;                                        //
    private:
        void                        rebuildQMatrix();                                       //
        double                      currentStationaryPrior = 0;                             //
        double                      kLambda;                                                //
        double                      oldStationaryPrior = 0;                                 //
        double                      omegaLambda;                                            //
        double                      rLambda;                                                //
        int                         currentActiveOmegas = 1;                                //
        int                         kAcceptCount = 0;                                       //
        int                         kCount = 0;                                             //
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
};

#endif