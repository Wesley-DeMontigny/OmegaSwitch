#ifndef EXPONENTIAL_PRIOR_HPP
#define EXPONENTIAL_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"


class ExponentialPrior : public PriorNode {
    public:
        ExponentialPrior(void)=delete;
        ExponentialPrior(BasicParameter<double>* l, BasicParameter<double>* p);
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();
        std::string writeValue() {return std::to_string(currentLnPrior);}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        BasicParameter<double>* lambda;
        BasicParameter<double>* param;
};

#endif