#ifndef PRIOR_NODE_HPP
#define PRIOR_NODE_HPP
#include "modeling/ModelNode.hpp"
#include <vector>
#include <string>

class PriorNode : public ModelNode {
    public:
        virtual void accept()=0;
        virtual void reject()=0;
        virtual void regenerate()=0;
        virtual std::string writeValue()=0;
        virtual double lnPrior()=0;
        virtual void sample()=0;
    protected:
        double currentLnPrior;
        double oldLnPrior;
};

#endif