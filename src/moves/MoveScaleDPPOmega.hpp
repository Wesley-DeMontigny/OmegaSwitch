#ifndef MOVE_SCALE_DPP_OMEGA_HPP
#define MOVE_SCALE_DPP_OMEGA_HPP
#include "Move.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"

class MoveScaleDPPOmega : public Move {
    public:
        MoveScaleDPPOmega(DirichletProcessPrior* d);
        double update();
        void tune();
    private:
        DirichletProcessPrior* dpp;
        double delta;
};

#endif