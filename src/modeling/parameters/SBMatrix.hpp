#ifndef SB_MATRIX_HPP
#define SB_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class SBMatrix : public Parameter {
    public:
        SBMatrix(Settings settings);
        Matrix<double> Q();
        std::vector<double> stationary();
        std::vector<double> dNdS();

        void accept();
        void reject();
        void tune();
        double lnPrior();

        double updateK();
        double updateOmega();
        double updateR();
        double updateStationary();
        double updateProportions();

        std::vector<double> getStationary();
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        std::vector<double> getOmegas() {return currentOmegas;}
        std::vector<double> getProportions() {return proportions;}
        double getR() {return currentR;}

        double kDelta;
        double omegaDelta;
        double rDelta;
        double stationaryAlpha;
        double proportionAlpha;
        int kCount = 0;
        int kAcceptCount = 0;
        int stationaryCount = 0;
        int stationaryAcceptCount = 0;
        int omegaCount = 0;
        int omegaAcceptCount = 0;
        int rCount = 0;
        int rAcceptCount = 0;
        int proportionCount = 0;
        int proportionAcceptCount = 0;
    private:
        int truncation;
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        void rebuildQMatrix();

        int moveChoice = -1;

        std::vector<int> randomStates;

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

        std::vector<double> currentOmegas;
        std::vector<double> oldOmegas;
        double currentOmegasPrior = 0;
        double oldOmegasPrior = 0;

        std::vector<double> proportions;
        std::vector<double> currentProportionParams;
        std::vector<double> oldProportionParams;
        //double currentProportionsPrior = 0;
        //double oldProportionsPrior = 0;

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