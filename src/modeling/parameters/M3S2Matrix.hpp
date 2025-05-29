#ifndef M3S2_MATRIX_HPP
#define M3S2_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class M3S2Matrix : public Parameter {
    public:
        M3S2Matrix(Settings settings);
        Matrix<double> Q();
        std::vector<double> stationary();
        void accept();
        void reject();
        double updateK();
        double updateOmega1();
        double updateOmega2();
        double updateOmega3();
        double updateGamma();
        double updateR1();
        double updateR2();
        double updateStationary();
        void tune();
        double lnPrior();
        std::vector<double> getStationary();
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getOmega1() {return currentOmega1;}
        double getOmega2() {return currentOmega2;}
        double getOmega3() {return currentOmega3;}
        double getGamma() {return currentGamma;}
        double getR1() {return currentR1;}
        double getR2() {return currentR2;}

        int kCount = 0;
        int kAcceptCount = 0;
        int stationaryCount = 0;
        int stationaryAcceptCount = 0;
        int omega1Count = 0;
        int omega1AcceptCount = 0;
        int omega2Count = 0;
        int omega2AcceptCount = 0;
        int omega3Count = 0;
        int omega3AcceptCount = 0;
        int gammaCount = 0;
        int gammaAcceptCount = 0;
        int r1Count = 0;
        int r1AcceptCount = 0;
        int r2Count = 0;
        int r2AcceptCount = 0;
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        void rebuildQMatrix();

        int moveChoice = -1;

        double kDelta;
        double omega1Delta;
        double omega2Delta;
        double omega3Delta;
        double gammaDelta;
        double r1Delta;
        double r2Delta;
        double stationaryAlpha;

        std::vector<int> randomStates;

        double kLambda;
        double gammaLambda;
        double rLambda;
        double omegaLambda;

        double currentK = 0;
        double oldK = 0;
        double currentKPrior = 0;
        double oldKPrior = 0;

        double currentGamma = 0;
        double oldGamma = 0;
        double currentGammaPrior = 0;
        double oldGammaPrior = 0;

        double currentR1 = 0;
        double oldR1 = 0;
        double currentR1Prior = 0;
        double oldR1Prior = 0;

        double currentR2 = 0;
        double oldR2 = 0;
        double currentR2Prior = 0;
        double oldR2Prior = 0;

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