#ifndef MOVE_TREE_NNI_HPP
#define MOVE_TREE_NNI_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "Move.hpp"

class MoveTreeNNI : public Move{
    public:
        MoveTreeNNI(TreeParameter* t);
        double update();
        void tune();
    private:
        TreeParameter* param;
};

#endif