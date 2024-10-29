#include "MoveScheduler.hpp"
#include "core/RandomVariable.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <iostream>

MoveScheduler::MoveScheduler(void) : totalWeight(0.0) {}

double MoveScheduler::updateRandom(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double pick = rng.uniformRv() * totalWeight;

    double cumulative = 0.0;

    for(int i = 0; i < params.size(); i++){
        if(pick - cumulative < weights[i]){
            return params[i]->update();
        }
        cumulative += weights[i];
    }

    Msg::error("The Move Handler did not find any parameter!");

    return 0;
}

void MoveScheduler::registerParam(Parameter* p, double weight){
    params.push_back(p);
    weights.push_back(weight);
    totalWeight += weight;
}