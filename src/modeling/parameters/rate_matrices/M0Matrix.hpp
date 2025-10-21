#ifndef M0_MATRIX_HPP
#define M0_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>
#include <array>

class Settings;

/**
 * @brief 
 * 
 */
class M0Matrix : public Parameter {
    public:
                                M0Matrix(Settings& settings);                   //

        double                  dNdS();                                         //
        double                  getK() {return currentParams[0];}               //
        double                  getKRate();                                     //
        double                  getOmega() {return currentParams[1];}           //
        double                  getOmegaRate();                                 //
        double                  getStationaryRate();                            //
        double                  lnPrior();                                      //
        double                  updateK();                                      //
        double                  updateOmega();                                  //
        double                  updateStationary();                             //
        Matrix<double>          Q();                                            //
        std::vector<double>     getRawStationary() {return currentStationary;}  //
        std::vector<double>     getStationary();                                //
        void                    accept();                                       //
        void                    reject();                                       //
        void                    tune();                                         //

        double                  kDelta;                                         //
        double                  omegaDelta;                                     //
        double                  stationaryAlpha;                                //
    private:
        void                    rebuildQMatrix();                               //

        double                  oldStationaryPrior = 0;                         //
        double                  currentStationaryPrior = 0;                     //
        double                  kLambda;                                        //
        double                  omegaLambda;                                    //
        int                     kAcceptCount = 0;                               //
        int                     kCount = 0;                                     //
        int                     moveChoice = -1;                                //
        int                     omegaAcceptCount = 0;                           //
        int                     omegaCount = 0;                                 //
        int                     stationaryAcceptCount = 0;                      //
        int                     stationaryCount = 0;                            //
        Matrix<double>          currentQMatrix;                                 //
        Matrix<double>          oldQMatrix;                                     //
        std::array<double, 2>   currentParamPriors = {0,0};                     //
        std::array<double, 2>   currentParams = {0,0};                          // Parameters for the rate matrix in the order of K, O
        std::array<double, 2>   oldParamPriors = {0,0};                         //
        std::array<double, 2>   oldParams = {0,0};                              //
        std::vector<double>     stationaryPriorAlpha;                           //
        std::vector<double>     currentStationary;                              //
        std::vector<double>     oldStationary;                                  //
        std::vector<int>        randomStates;                                   //
};

#endif