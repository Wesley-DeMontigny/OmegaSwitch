#include "MoveScaleDouble.hpp"
#include "core/RandomVariable.hpp"
#include <iostream>
#include <cmath>

MoveScaleDouble::MoveScaleDouble(BasicParameter<double>* d, double lower, double upper) : param({d}), delta(0.25), lowerBound(lower), upperBound(upper) {}

MoveScaleDouble::MoveScaleDouble(std::vector<BasicParameter<double>*> d, double lower, double upper) : param(d), delta(0.25), lowerBound(lower), upperBound(upper) {}
        
double MoveScaleDouble::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int randomMember = (int)(rng.uniformRv() * param.size());

    double newV = 0.0;
    double scale = 0.0;

    do {
        scale = std::exp(delta * (rng.uniformRv() - 0.5));
        newV = param[randomMember]->getValue() * scale;
    }
    while(newV <= lowerBound || newV >= upperBound);

    param[randomMember]->setValue(newV);
    param[randomMember]->dirty();

    return scale;
}

void MoveScaleDouble::tune(){
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