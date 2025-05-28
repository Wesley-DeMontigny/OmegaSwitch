#ifndef M1_MATRIX_HPP
#define M1_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class M1Matrix : public Parameter {
    public:
        M1Matrix(Settings settings);
        Matrix<double> Q();
        std::vector<double> stationary();
        void accept();
        void reject();
        double updateK();
        double updateOmega();
        double updateStationary();
        void tune();
        double lnPrior();
        double kPrior() {return currentKPrior;}
        double omegaPrior() {return currentOmega;}
        double stationaryPrior() {return currentStationaryPrior;}
        std::vector<double> getStationary();
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getOmega() {return currentOmega;}

        int kCount;
        int kAcceptCount;
        int stationaryCount;
        int stationaryAcceptCount;
        int omegaCount;
        int omegaAcceptCount;
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        int moveChoice;
        double kDelta;
        double omegaDelta;
        double stationaryAlpha;

        std::vector<int> randomStates;

        double kLambda;
        double currentK;
        double oldK;
        double currentKPrior;
        double oldKPrior;

        double omegaLambda;
        double currentOmega;
        double oldOmega;
        double currentOmegaPrior;
        double oldOmegaPrior;

        std::vector<double> currentStationary;
        std::vector<double> oldStationary;
        double currentStationaryPrior;
        double oldStationaryPrior;
        std::vector<double> stationaryPriorAlpha;

        std::set<std::pair<int, int>> nonsynonymous;
        std::set<std::pair<int, int>> synonymous;
        std::set<std::pair<int, int>> valid;
        std::set<std::pair<int, int>> transition;
};

#endif