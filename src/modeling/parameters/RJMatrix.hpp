#ifndef RJ_MATRIX_HPP
#define RJ_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class RJMatrix : public Parameter {
    public:
        RJMatrix(Settings settings);
        Matrix<double> Q();
        std::vector<double> stationary();
        std::tuple<double, double, double, double, double> dNdS();

        void accept();
        void reject();
        void tune();
        double lnPrior();

        double updateK();
        double updateOmega();
        double updateR();
        double updateStationary();
        double updateActiveOmegas();

        std::vector<double> getStationary();
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getOmega1() {return currentOmega1;}
        double getOmega2() {return currentOmega2;}
        double getOmega3() {return currentOmega3;}
        double getOmega4() {return currentOmega4;}
        double getOmega5() {return currentOmega5;}
        double getR() {return currentR;}
        double getActiveOmegas() {return currentActiveOmegas;}
        void setActiveOmegas(int o) {currentActiveOmegas = o; rebuildQMatrix();}

        double kDelta;
        double omegaDelta;
        double rDelta;
        double stationaryAlpha;

        int kCount = 0;
        int kAcceptCount = 0;
        int stationaryCount = 0;
        int stationaryAcceptCount = 0;
        int omegaCount = 0;
        int omegaAcceptCount = 0;
        int rCount = 0;
        int rAcceptCount = 0;
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        void rebuildQMatrix();

        int moveChoice = -1;

        std::vector<int> randomStates;

        int currentActiveOmegas = 5;
        int oldActiveOmegas = 5;

        double kLambda;
        double rLambda;
        double omegaLambda;

        double currentK = 0;
        double oldK = 0;
        double currentKPrior = 0;
        double oldKPrior = 0;

        double currentR = 0;
        double oldR = 0;
        double currentRPrior = 0;
        double oldRPrior = 0;

        double currentOmega1 = 0;
        double oldOmega1 = 0;
        double currentOmega1Prior = 0;
        double oldOmega1Prior = 0;

        double currentOmega2 = 0;
        double oldOmega2 = 0;
        double currentOmega2Prior = 0;
        double oldOmega2Prior = 0;

        double currentOmega3 = 0;
        double oldOmega3 = 0;
        double currentOmega3Prior = 0;
        double oldOmega3Prior = 0;

        double currentOmega4 = 0;
        double oldOmega4 = 0;
        double currentOmega4Prior = 0;
        double oldOmega4Prior = 0;

        double currentOmega5 = 0;
        double oldOmega5 = 0;
        double currentOmega5Prior = 0;
        double oldOmega5Prior = 0;

        std::vector<double> currentStationary;
        std::vector<double> oldStationary;
        double currentStationaryPrior = 0;
        double oldStationaryPrior = 0;
        std::vector<double> stationaryPriorAlpha;

        std::set<std::pair<int, int>> nonsynonymous;
        std::set<std::pair<int, int>> synonymous;
        std::set<std::pair<int, int>> valid;
        std::set<std::pair<int, int>> transition;
};

#endif