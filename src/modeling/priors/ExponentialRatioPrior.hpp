#ifndef EXPONENTIAL_RATIO_PRIOR_HPP
#define EXPONENTIAL_RATIO_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"

// This was a hard distribution to track down?
// See Field Guide to Continuous Probability Distributions by GE Crooks (pg. 40)
// This is apparently a form of the Pareto distribution
// Here I let s = 1
class ExponentialRatioPrior : public PriorNode {
    public:
        ExponentialRatioPrior(void)=delete;
        ExponentialRatioPrior(BasicParameter<double>* p);
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();
        std::string writeValue() {return std::to_string(currentLnPrior);}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        double lnPdf(double v);
        BasicParameter<double>* param;
};

#endif