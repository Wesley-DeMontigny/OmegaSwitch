#ifndef MOVE_SCALE_BRANCH_HPP
#define MOVE_SCALE_BRANCH_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "Move.hpp"

class MoveScaleBranch : public Move{
    public:
        MoveScaleBranch(TreeParameter* t);
        double update();
        void tune();
    private:
        TreeParameter* param;
        double delta;
};

#endif