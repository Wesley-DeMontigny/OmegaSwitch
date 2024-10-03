#include "MoveScaleDPPCategory.hpp"
#include "core/RandomVariable.hpp"
#include <iostream>
#include <cmath>

MoveScaleDPPCategory::MoveScaleDPPCategory(DirichletProcessPrior* d) : dpp(d), delta(0.25) {}

double MoveScaleDPPCategory::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int randomCategory = (int)(rng.uniformRv() * dpp->getNumCategories());

    double scale = std::exp(delta * (rng.uniformRv() - 0.5));
    double newV = dpp->getCategoryValue(randomCategory) * scale;

    dpp->setCategoryValue(randomCategory, newV);

    dpp->dirty();

    return scale;
}

void MoveScaleDPPCategory::tune(){
    double rate = (double)acceptedSinceTune/(double)countSinceTune;

    if ( rate > 0.44 ) {
        delta *= (1.0 + ((rate-0.44)/0.766));
    }
    else {
        delta /= (2.0 - rate/0.44);
    }

    acceptedSinceTune = 0;
    countSinceTune = 0;
}