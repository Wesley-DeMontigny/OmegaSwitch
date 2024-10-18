#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <string>

class TreeParameter : public Parameter{
    public:
        TreeParameter(void)=delete;
        TreeParameter(Alignment* aln, double lambda);
        ~TreeParameter();
        TreeParameter& operator=(const TreeParameter& t);
        TreeObject* getTree(){return trees[0];}
        const TreeObject* getTreeConst() const {return trees[0];}
        std::vector<double> getBranchLengths();

        void accept();
        void reject();

        double update();
        double lnPrior();

        std::string writeNewick() {trees[0]->getNewick();}
    private:
        double lambda;
        double currentPrior;
        double oldPrior;
        TreeObject* trees[2];
};

#endif