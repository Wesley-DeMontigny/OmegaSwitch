#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <string>

/**
 * @brief 
 * 
 */
class TreeParameter : public Parameter{
    public:
                            TreeParameter(void)=delete;
                            TreeParameter(Alignment* aln, std::string newick, double lambda);   //
                            TreeParameter(TreeObject& tree, double lambda);                     //
                            ~TreeParameter();                                                   //

        double              getBranchRate();                                                    //
        double              getTreeRate();                                                      //
        double              lnPrior();                                                          //
        double              update();                                                           //
        std::string         writeNewick() {return trees[0]->getNewick();}                       //
        TreeObject*         getTree(){return trees[0];}                                         //
        void                accept();                                                           //
        void                reject();                                                           //
        void                tune();                                                             //

        double              branchDelta;                                                        //
        double              treeAlpha;                                                          //
    private:
        bool                fixedTree;                                                          //
        double              currentPrior;                                                       //
        double              lambda;                                                             //
        double              oldPrior;                                                           //
        int                 branchAcceptCount;                                                  //
        int                 branchCount;                                                        //
        int                 moveChoice;                                                         //
        int                 treeAcceptCount;                                                    //
        int                 treeCount;                                                          //
        TreeObject*          trees[2];                                                          //
};

#endif