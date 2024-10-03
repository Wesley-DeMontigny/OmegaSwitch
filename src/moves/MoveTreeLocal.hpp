#ifndef MOVE_TREE_LOCAL_HPP
#define MOVE_TREE_LOCAL_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "Move.hpp"
#include <cmath>

class MoveTreeLocal : public Move{
    public:
        MoveTreeLocal(TreeParameter* t);
        double update();
        void tune();
    private:
        TreeParameter* param;
        double delta = std::log(4.0);
};

#endif