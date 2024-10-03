#include "MoveScheduler.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/RandomVariable.hpp"
#include "core/Msg.hpp"
#include "Move.hpp"
#include <iostream>

MoveScheduler::MoveScheduler(void) : totalWeight(0.0) {}

MoveScheduler::MoveScheduler(std::vector<Move*> moves, std::vector<double> w) : moveHandlers(moves), weights(w), totalWeight(0.0) {
    for(double i : weights)
        totalWeight += i;
}

// Just pick a random move for now 
Move* MoveScheduler::getMove(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double pick = rng.uniformRv() * totalWeight;

    double cumulative = 0.0;

    for(int i = 0; i < moveHandlers.size(); i++){
        if(pick - cumulative < weights[i]){
            return moveHandlers[i];
        }
        cumulative += weights[i];
    }

    Msg::error("The Move Handler did not find any move!");

    return nullptr;
}

void MoveScheduler::registerMove(Move* m, double weight){
    moveHandlers.push_back(m);
    weights.push_back(weight);
    totalWeight += weight;
}

void MoveScheduler::clearRecord(){
    for(Move* m : moveHandlers)
        m->clearRecord();
}

void MoveScheduler::tune(){
    for(Move* m : moveHandlers)
        m->tune();

}