#ifndef M0_MATRIX_HPP
#define M0_MATRIX_HPP
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
 * 
 */
class M0Matrix : public Parameter {
    public:
                                    M0Matrix(Settings& settings);                           //

        double                      dNdS() const;                                           //
        double                      getK() const {return currentParams[0];}                 //
        double                      getKRate() const;                                       //
        double                      getOmega() const {return currentParams[1];}             //
        double                      getOmegaRate() const;                                   //
        double                      getStationaryRate() const;                              //
        double                      lnPrior() override;                                     //
        double                      updateK();                                              //
        double                      updateOmega();                                          //
        double                      updateStationary();                                     //
        Matrix<double>              Q() const;                                              //
        const std::vector<double>&  getRawStationary() const {return currentStationary;}    //
        std::vector<double>         getStationary();                                        //
        void                        accept() override;                                      //
        void                        reject() override;                                      //
        void                        tune() override;                                        //

        double                      kDelta;                                                 //
        double                      omegaDelta;                                             //
        double                      stationaryAlpha;                                        //
    private:
        void                        rebuildQMatrix();                                       //
        double                      currentStationaryPrior = 0;                             //
        const double                kLambda;                                                //
        double                      oldStationaryPrior = 0;                                 //
        const double                omegaLambda;                                            //
        int                         kAcceptCount = 0;                                       //
        int                         kCount = 0;                                             //
        int                         omegaAcceptCount = 0;                                   //
        int                         omegaCount = 0;                                         //
        int                         stationaryAcceptCount = 0;                              //
        int                         stationaryCount = 0;                                    //
        Matrix<double>              currentQMatrix;                                         //
        Matrix<double>              oldQMatrix;                                             //
        MatrixMoves                 moveChoice = MatrixMoves::NO_MOVE;                      //
        std::array<double, 2>       currentParamPriors = {0,0};                             //
        std::array<double, 2>       currentParams = {0,0};                                  // Parameters for the rate matrix in the order of K, O
        std::array<double, 2>       oldParamPriors = {0,0};                                 //
        std::array<double, 2>       oldParams = {0,0};                                      //
        std::vector<double>         currentStationary;                                      //
        std::vector<double>         oldStationary;                                          //
        std::vector<double>         stationaryPriorAlpha;                                   //
        std::vector<int>            randomStates;                                           //
};

#endif