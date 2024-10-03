#ifndef MOVE_SCALE_DPP_HPP
#define MOVE_SCALE_DPP_HPP
#include "Move.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"

class MoveScaleDPPCategory : public Move {
    public:
        MoveScaleDPPCategory(DirichletProcessPrior* d);
        double update();
        void tune();
    private:
        DirichletProcessPrior* dpp;
        double delta;
};

#endif