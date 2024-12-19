#ifndef CODON_MULTI_MATRIX_HPP
#define CODON_MULTI_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class CodonMultiMatrix : public Parameter {
    public:
        CodonMultiMatrix(Settings settings, std::vector<double> pi);
        Matrix<double> Q(double omega1, double omega2);
        std::vector<double> stationary();
        void accept();
        void reject();
        double update();
        void tune();
        double lnPrior();
        double kPrior() {return currentKPrior;}
        double rPrior() {return currentRPrior;}
        bool updatingStationary() {return updateStationary;}
        std::vector<double> getStationary();
        std::vector<double> getRawStationary() {return currentStationary;}
        double getK() {return currentK;}
        double getR() {return currentR;}
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        int moveChoice;
        int kCount;
        int kAcceptCount;
        double kDelta;
        int rCount;
        int rAcceptCount;
        double rDelta;
        int stationaryCount;
        int stationaryAcceptCount;
        double stationaryAlpha;
        double stationaryOffset;

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
        bool updateStationary;

        std::set<std::
        pair<int, int>> nonsynonymous;
        std::set<std::pair<int, int>> valid;
        std::set<std::pair<int, int>> transition;

        std::vector<int> randomStates;
};

#endif