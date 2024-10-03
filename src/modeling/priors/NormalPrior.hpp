#ifndef NORMAL_PRIOR_HPP
#define NORMAL_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"

class NormalPrior : public PriorNode {
    public:
        NormalPrior(void)=delete;
        NormalPrior(BasicParameter<double>* m, BasicParameter<double>* s, BasicParameter<double>* p);
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();
        std::string writeValue() {return std::to_string(currentLnPrior);}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        BasicParameter<double>* mu;
        BasicParameter<double>* sigma;
        BasicParameter<double>* param;
};

#endif