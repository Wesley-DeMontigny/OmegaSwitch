#include "MoveDPPBeta.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include <iostream>
#include <cmath>

MoveDPPBeta::MoveDPPBeta(DirichletProcessPrior* d) : dpp(d), alpha(1.0) {}

double MoveDPPBeta::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    int randomCategory = (int)(rng.uniformRv() * dpp->getNumCategories());
    double betaVal = dpp->getCategoryBeta(randomCategory);

    double a = alpha + 1.0;
    double b = (alpha / betaVal) - a + 2.0;
    double newVal = Probability::Beta::rv(&rng, a, b);

    dpp->setCategoryBeta(randomCategory, newVal);

    double scalingFactor = (1.0 - newVal)/(1.0 - betaVal);

    double forward = Probability::Beta::lnPdf(a, b, newVal);
    double newA = alpha + 1.0;
    double newB = (alpha / newVal) - a + 2.0;
    double backward = Probability::Beta::lnPdf(newA, newB, betaVal);
    
    double hastings = backward - forward;

    dpp->dirty();

    return hastings;
}

void MoveDPPBeta::tune(){
    double rate = (double)acceptedSinceTune/(double)countSinceTune;

    if ( rate > 0.44 ) {
        alpha *= (1.0 + ((rate-0.44)/0.766));
    }
    else {
        alpha /= (2.0 - rate/0.44);
    }

    acceptedSinceTune = 0;
    countSinceTune = 0;
}