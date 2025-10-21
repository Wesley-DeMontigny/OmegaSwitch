#ifndef RJ_MATRIX_HPP
#define RJ_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>
#include <array>

class Settings;

/**
 * @brief 
 */
class RJMatrix : public Parameter {
    public:
                                RJMatrix(Settings& settings);                           //

        double                  getActiveOmegas() {return currentActiveOmegas;}         //
        double                  getK() {return currentParams[0];}                       //
        double                  getKRate();                                             //
        double                  getOmega(int i) {return currentParams[2+i];}            // Get the ith Omega (or Omega increment)
        double                  getOmegaRate();                                         //
        double                  getR() {return currentParams[1];}                       //
        double                  getRRate();                                             //
        double                  getStationaryRate();                                    //
        double                  lnPrior();                                              //
        double                  updateActiveOmegas();                                   //
        double                  updateK();                                              //
        double                  updateOmega();                                          //
        double                  updateR();                                              //
        double                  updateStationary();                                     //
        Matrix<double>          Q();                                                    //
        std::array<double, 5>   dNdS();                                                 //
        std::vector<double>     getRawStationary() {return currentStationary;}          //
        std::vector<double>     getStationary();                                        //
        void                    accept();                                               //
        void                    reject();                                               //
        void                    setActiveOmegas(int o) {                                //
                                currentActiveOmegas = o; 
                                rebuildQMatrix();
                                }
        void                    tune();                                                 //

        double                  kDelta;                                                 //
        double                  omegaDelta;                                             //
        double                  rDelta;                                                 //
        double                  stationaryAlpha;                                        //
    private:
        void                    rebuildQMatrix();                                       //

        double                  currentStationaryPrior = 0;                             //
        double                  kLambda;                                                //
        double                  oldStationaryPrior = 0;                                 //
        double                  omegaLambda;                                            //
        double                  rLambda;                                                //
        int                     currentActiveOmegas = 1;                                //
        int                     kAcceptCount = 0;                                       //
        int                     kCount = 0;                                             //
        int                     moveChoice = -1;                                        //
        int                     oldActiveOmegas = 1;                                    //
        int                     omegaAcceptCount = 0;                                   //
        int                     omegaCount = 0;                                         //
        int                     rAcceptCount = 0;                                       //
        int                     rCount = 0;                                             //
        int                     stationaryAcceptCount = 0;                              //
        int                     stationaryCount = 0;                                    //
        Matrix<double>          currentQMatrix;                                         //
        Matrix<double>          oldQMatrix;                                             //
        std::array<double, 7>   currentParamPriors = {0,0,0,0,0,0,0};                   //
        std::array<double, 7>   currentParams = {0,0,0,0,0,0,0};                        // Parameters for the rate matrix in the order of K, R, O1, O2, O3, O4, O5
        std::array<double, 7>   oldParamPriors = {0,0,0,0,0,0,0};                       //
        std::array<double, 7>   oldParams = {0,0,0,0,0,0,0};                            // 
        std::vector<double>     currentStationary;                                      // 
        std::vector<double>     oldStationary;                                          //  
        std::vector<double>     stationaryPriorAlpha;                                   // 
        std::vector<int>        randomStates;                                           // 
};

#endif