#define _USE_MATH_DEFINES
#include "NormalPrior.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include <cmath>
#include <iostream>

NormalPrior::NormalPrior(BasicParameter<double>* m, BasicParameter<double>* s, BasicParameter<double>* p) : currentLnPrior(0.0), oldLnPrior(0.0), mu(m), sigma(s), param(p) {
    this->dirty();
}

void NormalPrior::accept() {
    oldLnPrior = currentLnPrior;

    mu->accept();
    mu->clean();
    sigma->accept();
    sigma->clean();
    param->accept();
    param->clean();
}

void NormalPrior::reject() {
    currentLnPrior = oldLnPrior;

    mu->reject();
    mu->clean();
    sigma->reject();
    sigma->clean();
    param->reject();
    param->clean();
}

void NormalPrior::regenerate(){
    mu->regenerate();
    sigma->regenerate();
    param->regenerate();

    if(mu->isDirty() || sigma->isDirty() || param->isDirty())
        this->dirty();
    

    if(this->isDirty()){
        currentLnPrior = Probability::Normal::lnPdf(mu->getValue(), sigma->getValue(), param->getValue());
    }
}

void NormalPrior::sample(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    param->setValue(Probability::Normal::rv(&rng, mu->getValue(), sigma->getValue()));
    param->accept();
    param->regenerate();
}