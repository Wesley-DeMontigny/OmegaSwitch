#ifndef GAMMA_PRIOR_HPP
#define GAMMA_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"
#include <vector>

class GammaPrior : public PriorNode {
    public:
        GammaPrior(void)=delete;
        GammaPrior(BasicParameter<double>* a, BasicParameter<double>* b, BasicParameter<double>* p);
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();
        std::string writeValue() {return std::to_string(currentLnPrior);}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        BasicParameter<double>* shape;
        BasicParameter<double>* rate;
        BasicParameter<double>* param;
};

#endif