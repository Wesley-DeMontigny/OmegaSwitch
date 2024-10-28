#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <string>

class TreeParameter : public Parameter{
    public:
        TreeParameter(void)=delete;
        TreeParameter(Alignment* aln, std::string newick, double lambda);
        ~TreeParameter();
        TreeObject* getTree(){return trees[0];}

        void accept();
        void reject();

        double update();
        void tune();
        double lnPrior();

        std::string writeNewick() {return trees[0]->getNewick();}
    private:
        bool fixedTree;
        int moveChoice;
        int count;
        int acceptCount;
        double delta;
        double lambda;
        double currentPrior;
        double oldPrior;
        TreeObject* trees[2];
};

#endif