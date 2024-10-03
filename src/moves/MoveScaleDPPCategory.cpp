#include "MoveScaleDPPCategory.hpp"
#include "core/RandomVariable.hpp"
#include <iostream>
#include <cmath>

MoveScaleDPPCategory::MoveScaleDPPCategory(DirichletProcessPrior* d) : dpp(d), delta(0.25) {}

double MoveScaleDPPCategory::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int randomCategory = (int)(rng.uniformRv() * dpp->getNumCategories());

    double scale = std::exp(delta * (rng.uniformRv() - 0.5));
    
    bool omega1 = rng.uniformRv() > 0.5;

    if(omega1) {
        double newV = dpp->getCategoryOmega1(randomCategory) * scale;
        dpp->setCategoryOmega1(randomCategory, newV);
    }
    else {
        double newV = dpp->getCategoryOmega2(randomCategory) * scale;
        dpp->setCategoryOmega2(randomCategory, newV);
    }

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