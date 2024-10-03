#include "MoveBetaSimplex.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include <iostream>
#include <cmath>

MoveBetaSimplex::MoveBetaSimplex(std::vector<BasicParameter<double>*> d) : params(d), alpha(1.0) {}
        
double MoveBetaSimplex::update(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int paramNum = params.size();
    int i = (int)(rng.uniformRv() * paramNum);

    BasicParameter<double>* p = params[i];
    double pVal = p->getValue();

    double a = alpha + 1.0;
    double b = (alpha / pVal) - a + 2.0;
    double newVal = Probability::Beta::rv(&rng, a, b);

    p->setValue(newVal);

    double scalingFactor = (1.0 - newVal)/(1.0 - pVal);

    double sum = 0.0;
    double hastings = 0.0;
    for(BasicParameter<double>* param : params){
        if(param != p)
            param->setValue(param->getValue() * scalingFactor);

        if(param->getValue() <= 1E-100)
            return -1.0 * INFINITY;
        
        sum += param->getValue();
    }

    //Normalize to make sure this doesn't drift from 1.0;
    for (BasicParameter<double>* param : params) {
        double paramVal = param->getValue();
        param->setValue(paramVal/sum);

        param->dirty();
    }

    // The probability of getting our new value
    double forward = Probability::Beta::lnPdf(a, b, newVal);
    double newA = alpha + 1.0;
    double newB = (alpha / newVal) - a + 2.0;
    // The probability of getting our old value in the future
    double backward = Probability::Beta::lnPdf(newA, newB, pVal);
    
    hastings = backward - forward;
    
    hastings += (paramNum - 2) * std::log(scalingFactor) - (paramNum - 1) * std::log(sum);

    return hastings;
}

void MoveBetaSimplex::tune(){
    double rate = (double)acceptedSinceTune/(double)countSinceTune;

    if (rate > 0.44) {
        alpha /= (1.0 + ((rate-0.44)/(1.0 - 0.44)));
    }
    else {
        alpha *= (2.0 - rate/0.44);
    }

    alpha = fmin(alpha, 100);

    acceptedSinceTune = 0;
    countSinceTune = 0;
}