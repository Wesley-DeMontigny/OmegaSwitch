#ifndef MOVE_SCHEDULER_HPP
#define MOVE_SCHEDULER_HPP
#include <vector>

class Parameter;

//For now this is hardset to only be moves on the tree
class MoveScheduler {
    public:
        MoveScheduler(void);
        double updateRandom();
        void registerParam(Parameter* m, double weight);
    private:
        std::vector<Parameter*> params;
        double totalWeight;
        std::vector<double> weights;

};

#endif