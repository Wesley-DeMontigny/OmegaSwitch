#include "DirichletPrior.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Msg.hpp"
#include <cmath>
#include <iostream>

DirichletPrior::DirichletPrior(std::vector<double> d, std::vector<BasicParameter<double>*> p) : currentLnPrior(0.0), oldLnPrior(0.0), dist(d), param(p) {
    if(dist.size() != param.size())
        Msg::error("The number of parameters must match the number of elements in the distribution!");

    this->dirty();
}

void DirichletPrior::accept() {
    oldLnPrior = currentLnPrior;

    for(BasicParameter<double>* p : param){
        p->accept();
        p->clean();
    }
}

void DirichletPrior::reject() {
    currentLnPrior = oldLnPrior;

    for(BasicParameter<double>* p : param){
        p->reject();
        p->clean();
    }  
}

void DirichletPrior::regenerate(){
    for(BasicParameter<double>* p : param){
        p->regenerate();
        if(p->isDirty())
            this->dirty();
    }
    

    if(this->isDirty()){
        std::vector<double> values;
        
        for(int i = 0; i < param.size(); i++){
            double val = param[i]->getValue();
            values.push_back(val);
        }
        currentLnPrior = Probability::Dirichlet::lnPdf(dist, values);
    }
}

void DirichletPrior::sample(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<double> randomVal = dist;

    Probability::Dirichlet::rv(&rng, dist, randomVal);
    for(int i = 0; i < param.size(); i++){
        param[i]->setValue(randomVal[i]);
        param[i]->accept();
        param[i]->regenerate();
    }
}