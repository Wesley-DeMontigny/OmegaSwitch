#ifndef MOVE_DPP_BETA_HPP
#define MOVE_DPP_BETA_HPP
#include "Move.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"

class MoveDPPBeta : public Move {
    public:
        MoveDPPBeta(DirichletProcessPrior* d);
        double update();
        void tune();
    private:
        DirichletProcessPrior* dpp;
        double alpha;
};

#endif