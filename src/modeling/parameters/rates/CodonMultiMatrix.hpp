#ifndef CODON_MULTI_MATRIX_HPP
#define CODON_MULTI_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"
#include "modeling/parameters/BasicParameter.hpp"
#include <set>

class CodonMultiMatrix : public ModelNode {
    public:
        CodonMultiMatrix(DirichletProcessPrior* d, BasicParameter<double>* k, BasicParameter<double>* r, std::vector<BasicParameter<double>*> pi);
        Matrix<double> Q(double omega1, double omega2);
        std::vector<double> stationary();
        void accept();
        void reject();
        void regenerate();
        std::string writeValue() {return "";}
    private:

        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        BasicParameter<double>* kParam;
        BasicParameter<double>* rParam;

        std::set<std::pair<int, int>> nonsynonymous;
        std::set<std::pair<int, int>> valid;
        std::set<std::pair<int, int>> transition;

        std::vector<BasicParameter<double>*> stationaryDist;
        std::vector<double> returnStationary;
};

#endif