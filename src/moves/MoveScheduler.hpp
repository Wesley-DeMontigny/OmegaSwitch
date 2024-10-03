#ifndef MOVE_SCHEDULER_HPP
#define MOVE_SCHEDULER_HPP
#include <vector>

class TreeParameter;
class Move;

//For now this is hardset to only be moves on the tree
class MoveScheduler {
    public:
        MoveScheduler(void);
        MoveScheduler(std::vector<Move*> moves, std::vector<double> w);
        Move* getMove();
        void registerMove(Move* m, double weight);
        //void registerMove(AbstractMove* m, int w);
        void tune();
        void clearRecord();
    private:
        std::vector<Move*> moveHandlers;
        double totalWeight;
        std::vector<double> weights;

};

#endif