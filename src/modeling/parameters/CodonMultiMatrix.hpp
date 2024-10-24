#ifndef CODON_MULTI_MATRIX_HPP
#define CODON_MULTI_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class CodonMultiMatrix : public Parameter {
    public:
        CodonMultiMatrix(double rL, double kL, std::vector<double> pi, bool updatePi=false);
        Matrix<double> Q(double omega1, double omega2);
        std::vector<double> stationary();
        void accept();
        void reject();
        double update();
        void tune();
        double lnPrior();
        double kPrior() {return currentKPrior;}
        double rPrior() {return currentRPrior;}
        double stationaryPrior() {return currentStationaryPrior;}
        bool updatingStationary() {return updateStationary;}
        std::vector<double> getStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getR() {return currentR;}
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        int moveChoice;
        int rCount;
        int rAcceptCount;
        double rDelta;
        int kCount;
        int kAcceptCount;
        double kDelta;
        int stationaryCount;
        int stationaryAcceptCount;
        double stationaryAlpha;

        double kLambda;
        double currentK;
        double oldK;
        double currentKPrior;
        double oldKPrior;

        double rLambda;
        double currentR;
        double oldR;
        double currentRPrior;
        double oldRPrior;

        std::vector<double> currentStationary;
        std::vector<double> oldStationary;
        double currentStationaryPrior;
        double oldStationaryPrior;
        bool updateStationary;
        std::vector<double> flatDirichlet;

        std::set<std::pair<int, int>> nonsynonymous;
        std::set<std::pair<int, int>> valid;
        std::set<std::pair<int, int>> transition;
};

#endif