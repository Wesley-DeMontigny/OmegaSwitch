#ifndef CODON_MULTI_MATRIX_HPP
#define CODON_MULTI_MATRIX_HPP
#include "RateMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/BasicParameter.hpp"
#include <set>

class DirichletProcessPrior;

class CodonMultiMatrix : public RateMatrix {
    public:
        CodonMultiMatrix(DirichletProcessPrior* o, BasicParameter<double>* k, std::vector<BasicParameter<double>*> pi);
        Matrix<double> Q() {return Matrix<double>(61, 61, 0.0);}
        Matrix<double> Q(int index);
        std::vector<double> stationary();
        void accept();
        void reject();
        void regenerate();
        std::string writeValue() {return "";}
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;
        DirichletProcessPrior* dNdS;
        BasicParameter<double>* transitionTransversionRatio;
        std::set<std::pair<int, int>> nonsynonymous;
        std::set<std::pair<int, int>> valid;
        std::set<std::pair<int, int>> transition;
        std::vector<BasicParameter<double>*> stationaryDist;
        std::vector<double> returnStationary;
};

#endif