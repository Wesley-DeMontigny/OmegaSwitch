#ifndef M0_MATRIX_HPP
#define M0_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class M0Matrix : public Parameter {
    public:
        M0Matrix(Settings settings);
        Matrix<double> Q();
        std::vector<double> stationary();
        void accept();
        void reject();
        double updateK();
        double updateOmega();
        double updateStationary();
        void tune();
        double lnPrior();
        std::vector<double> getStationary();
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getOmega() {return currentOmega;}

        int kCount = 0;
        int kAcceptCount = 0;
        int stationaryCount = 0;
        int stationaryAcceptCount = 0;
        int omegaCount = 0;
        int omegaAcceptCount = 0;
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        void rebuildQMatrix();

        int moveChoice = -1;
        double kDelta;
        double omegaDelta;
        double stationaryAlpha;

        std::vector<int> randomStates;

        double kLambda;
        double omegaLambda;

        double currentK = 0;
        double oldK = 0;
        double currentKPrior = 0;
        double oldKPrior = 0;

        double currentOmega = 0;
        double oldOmega = 0;
        double currentOmegaPrior = 0;
        double oldOmegaPrior = 0;

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