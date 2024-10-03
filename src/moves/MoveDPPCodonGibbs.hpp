#ifndef MOVE_DPP_MATRIX_GIBBS_HPP
#define MOVE_DPP_MATRIX_GIBBS_HPP
#include "Move.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"
#include "modeling/likelihoods/CodonMMPhyloCTMC.hpp"

class MoveDPPCodonGibbs : public Move {
    public:
        MoveDPPCodonGibbs(void)=delete;
        MoveDPPCodonGibbs(CodonMMPhyloCTMC* l, DirichletProcessPrior* d);
        double update();
        void tune();
    private:
        CodonMMPhyloCTMC* likelihood;
        DirichletProcessPrior* dpp;
        int currentMember;

};

#endif