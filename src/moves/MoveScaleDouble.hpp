#ifndef MOVE_SCALE_DOUBLE_HPP
#define MOVE_SCALE_DOUBLE_HPP
#include "Move.hpp"
#include "modeling/parameters/BasicParameter.hpp"

class MoveScaleDouble : public Move {
    public:
        MoveScaleDouble(BasicParameter<double>* d, double lower = -1 * INFINITY, double upper = INFINITY);
        MoveScaleDouble(std::vector<BasicParameter<double>*> d, double lower = -1 * INFINITY, double upper = INFINITY);
        double update();
        void tune();
    private:
        std::vector<BasicParameter<double>*> param;
        double delta;
        double lowerBound;
        double upperBound;
};

#endif