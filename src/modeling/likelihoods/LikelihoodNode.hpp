#ifndef LIKELIHOOD_NODE_HPP
#define LIKELIHOOD_NODE_HPP
#include "modeling/ModelNode.hpp"
#include <vector>
#include <string>

class LikelihoodNode : public ModelNode {
    public:
        virtual void accept()=0;
        virtual void reject()=0;
        virtual void regenerate()=0;
        virtual std::string writeValue()=0;
        virtual double lnLikelihood()=0;
    protected:
        double currentLnLikelihood;
        double oldLnLikelihood;
};

#endif