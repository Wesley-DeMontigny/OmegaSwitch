#ifndef DIRICHLET_PRIOR_HPP
#define DIRICHLET_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"


class DirichletPrior : public PriorNode {
    public:
        DirichletPrior(void)=delete;
        DirichletPrior(std::vector<double> d, std::vector<BasicParameter<double>*> p);
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();
        std::string writeValue() {return std::to_string(currentLnPrior);}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        std::vector<double> dist;
        std::vector<BasicParameter<double>*> param;
};

#endif