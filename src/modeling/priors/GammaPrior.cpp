#define _USE_MATH_DEFINES
#include "GammaPrior.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include <cmath>
#include <iostream>

GammaPrior::GammaPrior(BasicParameter<double>* a, BasicParameter<double>* b, BasicParameter<double>* p) : currentLnPrior(0.0), oldLnPrior(0.0), shape(a), rate(b), param(p) {
    this->dirty();
}

void GammaPrior::accept() {
    oldLnPrior = currentLnPrior;

    shape->accept();
    shape->clean();
    rate->accept();
    rate->clean();
    param->accept();
    param->clean();
}

void GammaPrior::reject() {
    currentLnPrior = oldLnPrior;

    shape->reject();
    shape->clean();
    rate->reject();
    rate->clean();
    param->reject();
    param->clean();
}

void GammaPrior::regenerate(){
    shape->regenerate();
    rate->regenerate();
    param->regenerate();

    if(shape->isDirty() || rate->isDirty() || param->isDirty())
        this->dirty();
    

    if(this->isDirty()){
        currentLnPrior = Probability::Gamma::lnPdf(shape->getValue(), rate->getValue(), param->getValue());
    }
}

void GammaPrior::sample(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    param->setValue(Probability::Gamma::rv(&rng, shape->getValue(), rate->getValue()));

    param->accept();
    param->regenerate();
}