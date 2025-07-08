#ifndef RJ_DPP_MATRIX_HPP
#define RJ_DPP_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class RJDPPMatrix : public Parameter {
    public:
        RJDPPMatrix(Settings settings);
        Matrix<double> Q(double omega1, double omega2, double omega3);
        std::vector<double> stationary();
        std::tuple<double, double, double> dNdS(double omega1, double omega2, double omega3);

        void accept();
        void reject();
        void tune();
        double lnPrior();

        double updateK();
        double updateR();
        double updateStationary();
        void refreshQBackground(int numClasses); // To be called to incorporate the stationary distribution, K, and R

        std::vector<double> getStationary(int omegaCount);
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getR() {return currentR;}

        double kDelta;
        double rDelta;
        double stationaryAlpha;

        int kCount = 0;
        int kAcceptCount = 0;
        int stationaryCount = 0;
        int stationaryAcceptCount = 0;
        int rCount = 0;
        int rAcceptCount = 0;
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        void rebuildQMatrix();

        int moveChoice = -1;

        std::vector<int> randomStates;

        double kLambda;
        double rLambda;

        double currentK = 0;
        double oldK = 0;
        double currentKPrior = 0;
        double oldKPrior = 0;

        double currentR = 0;
        double oldR = 0;
        double currentRPrior = 0;
        double oldRPrior = 0;

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