#ifndef RJ_DPP_MATRIX_HPP
#define RJ_DPP_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "core/Probability.hpp"
#include <set>
#include <vector>
#include <array>
#include <optional>

class Settings;

struct Category {
    std::array<double, 3> omegas;
    std::vector<int> members;
    bool dirty;
};

/**
 * @brief 
 * 
 */
class RJDPPMatrix : public Parameter {
    public:
                                RJDPPMatrix(Settings& settings, int nC);                //

        double                  getK() {return currentParams[0];}                       //
        double                  getKRate();                                             //
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
        double                  getAlpha() {return dpAlpha;}                            //
        double                  getActiveOmegas() {return currentActiveOmegas;}         //
        std::optional<int>      unassignMember(int member);                             //
        int                     getNumCategories() {return currentCategories.size();}   //
        Matrix<double>          Q(int i);                                               //
        Matrix<double>          Q(const std::array<double, 3>& omegas);                 //
        std::array<double, 3>   dNdS(int i);                                            //
        std::vector<Category>   getCategories() {return currentCategories;}             //
        std::vector<double>     getRawStationary() {return currentStationary;}          //
        std::vector<double>     getStationary(int omegaCount);                          //
        std::vector<int>        getAssignments() { return assignments;}                 //
        void                    accept();                                               //
        void                    reject();                                               //
        void                    setActiveOmegas(int o) {                                //
                                currentActiveOmegas = o; 
                                rebuildQMatrix();
                                }
        void                    tune();                                                 //
        void                    addCategory(const std::array<double, 3>& omegas);       //
        void                    assignMember(int member, int category);                 //
        void                    popCategories(int n);                                   //
        void                    regenerateAssignments();                                //
        void                    removeCategory(int index);                              //
        void                    regenerateCatPrior();                                   //

        double                  kDelta;                                                 // 
        double                  rDelta;                                                 //
        double                  omegaDelta;                                             //
        double                  stationaryAlpha;                                        //
    private:
        void                    rebuildQMatrix();                                       //
        double                  calculateAlpha(double expectedCat, int members);        //
        double                  expectedCategories(double a, int members);              //

        double                  currentStationaryPrior = 0;                             //
        double                  dpAlpha = 0.0;                                          //
        double                  kLambda;                                                //
        double                  oldStationaryPrior = 0;                                 //
        double                  omegaLambda;                                            //
        double                  rLambda;                                                //
        double                  currentCatPrior = 0.0;                                  //
        double                  oldCatPrior = 0.0;                                      //
        id_t                    omegaCount = 0;                                         //
        int                     currentActiveOmegas = 1;                                //
        int                     kAcceptCount = 0;                                       //
        int                     kCount = 0;                                             //
        int                     moveChoice = -1;                                        //
        int                     numChar;                                                //
        int                     oldActiveOmegas = 1;                                    //
        int                     omegaAcceptCount = 0;                                   //
        int                     rAcceptCount = 0;                                       //
        int                     rCount = 0;                                             //
        int                     stationaryAcceptCount = 0;                              //
        int                     stationaryCount = 0;                                    //
        Matrix<double>          currentQMatrix;                                         //
        Matrix<double>          oldQMatrix;                                             //
        std::array<double, 2>   currentParamPriors = {0,0};                             //
        std::array<double, 2>   currentParams = {0,0};                                  // Parameters for the rate matrix in the order of K, R
        std::array<double, 2>   oldParamPriors = {0,0};                                 //
        std::array<double, 2>   oldParams = {0,0};                                      //
        std::vector<Category>   currentCategories;                                      //
        std::vector<Category>   oldCategories;                                          //
        std::vector<double>     currentStationary;                                      //
        std::vector<double>     oldStationary;                                          //
        std::vector<double>     stationaryPriorAlpha;                                   //
        std::vector<int>        randomStates;                                           //
        std::vector<int>        assignments;                                            //
};

#endif