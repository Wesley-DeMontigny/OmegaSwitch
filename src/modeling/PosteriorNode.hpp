#ifndef POSTERIOR_NODE_HPP
#define POSTERIOR_NODE_HPP
#include "ModelNode.hpp"
#include "modeling/likelihoods/LikelihoodNode.hpp"
#include "modeling/priors/PriorNode.hpp"

class PosteriorNode : public ModelNode {
    public:
        PosteriorNode(void)=delete;
        PosteriorNode(LikelihoodNode* lN, std::vector<PriorNode*> pN);//May want to modify this to accept multiple priors and likelihoods
        void accept();
        void reject();
        void regenerate();
        double lnPosterior() {return currentLnPosterior;}
        std::string writeValue() {return std::to_string(currentLnPosterior);}
    private:
        LikelihoodNode* likelihood;
        std::vector<PriorNode*> priors;
        double currentLnPosterior;
        double oldLnPosterior;
};

#endif