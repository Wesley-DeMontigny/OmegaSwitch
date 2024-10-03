#ifndef MOVE_BETA_SIMPLEX_HPP
#define MOVE_BETA_SIMPLEX_HPP
#include "Move.hpp"
#include "modeling/parameters/BasicParameter.hpp"

class MoveBetaSimplex : public Move {
    public:
        MoveBetaSimplex(std::vector<BasicParameter<double>*> d);
        double update();
        void tune();
    private:
        std::vector<BasicParameter<double>*> params;
        double alpha;
};

#endif