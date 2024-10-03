#ifndef TREE_PRIOR_HPP
#define TREE_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"
#include <vector>
#include <set>

class TreeParameter;

class TreePrior : public PriorNode {
    public:
        TreePrior(void)=delete;
        TreePrior(TreeParameter* t);
        void setExponentialBranchPrior(BasicParameter<double>* l);
        void addMonophyleticConstraint(std::set<int> taxa, double strength);
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();
        std::string writeValue() {return std::to_string(currentLnPrior);}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        double numBranches;
        BasicParameter<double>* lambda;
        TreeParameter* tree;
        std::vector<std::set<int>> monophyleticContraints;
};

#endif