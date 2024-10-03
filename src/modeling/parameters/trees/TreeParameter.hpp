#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/ModelNode.hpp"
#include "TreeObject.hpp"
#include <string>

class TreeParameter : public ModelNode{
    public:
        TreeParameter(void)=delete;
        TreeParameter(int nt);
        TreeParameter(Alignment* aln);
        TreeParameter(std::string newick, std::vector<std::string> taxaNames);
        ~TreeParameter();
        TreeParameter& operator=(const TreeParameter& t);
        TreeObject* getTree(){return trees[0];}
        const TreeObject* getTreeConst() const {return trees[0];}
        std::vector<double> getBranchLengths();
        void accept();
        void reject();
        void regenerate();
        std::string writeValue();
    private:
        TreeObject* trees[2];
};

#endif