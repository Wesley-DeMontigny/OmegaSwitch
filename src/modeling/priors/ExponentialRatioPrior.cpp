#include "ExponentialRatioPrior.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include <cmath>
#include <iostream>

ExponentialRatioPrior::ExponentialRatioPrior(BasicParameter<double>* p) : currentLnPrior(0.0), oldLnPrior(0.0), param(p) {
    this->dirty();
}

void ExponentialRatioPrior::accept() {
    oldLnPrior = currentLnPrior;

    param->accept();
    param->clean();
}

void ExponentialRatioPrior::reject() {
    currentLnPrior = oldLnPrior;

    param->reject();
    param->clean();
}

void ExponentialRatioPrior::regenerate(){
    param->regenerate();

    if(param->isDirty())
            this->dirty();
    

    if(this->isDirty()){
        double val = param->getValue();
        if(val > 0)
            currentLnPrior = lnPdf(val);
        else
            currentLnPrior = -INFINITY;
    }
}

void ExponentialRatioPrior::sample(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    //Lambda value actually doesn't matter
    param->setValue(Probability::Exponential::rv(&rng, 1)/Probability::Exponential::rv(&rng, 1));
    param->accept();
    param->regenerate();
}

double ExponentialRatioPrior::lnPdf(double v){
    // log(1/(1+v)^2)
    return -2 * std::log(1.0 + v);
}