#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <string>

/**
 * @brief Enum for tracking the ID of the last proposed move
 */
enum TreeMoves{
    NO_MOVE = 1,
    BRANCH_MOVE = 2,
    TREE_MOVE = 3
};

/**
 * @brief The tree parameter object - as opposed to the tree object. The tree object is the 
 * underlying phylogenetic tree datatype, while the tree parameter type implements all the necessary
 * methods needed from a parameter that is being updated by MCMC. Thus, we need to contain two tree
 * objects so we can easily revert to the previous copy if an MCMC proposal has been reverted.
 */
class TreeParameter : public Parameter{
    public:
                            TreeParameter(void)=delete;
                            TreeParameter(Alignment& aln, std::string& newick, double lambda);  // Construct a new tree from an alignment, a newick string, and a rate parameter for the exponential prior
                            TreeParameter(TreeObject& tree, double lambda);                     // Construct a new tree from a pre-made tree object and a rate parameter for the exponential prior
                            ~TreeParameter();                                                   // Destructor

        double              getBranchRate() const;                                              // Get the acceptance rate for branch proposals
        double              getTreeRate() const;                                                // Get the acceptance rate for proposals on the whole tree lengths
        double              lnPrior() override;                                                 // Return the log prior probability of the tree
        double              update();                                                           // Propose a random MCMC move
        std::string         writeNewick() const {return trees[0]->getNewick();}                 // Return the newick string of the active tree object
        TreeObject*         getTree() const {return trees[0];}                                  // Get the active tree object
        void                accept() override;                                                  // Accept changes proposed to the tree
        void                reject() override;                                                  // Reject changes proposed to the tree
        void                tune() override;                                                    // Tune proposals that update the tree

        double              branchDelta;                                                        // The delta parameter for the scale branch move
        double              treeAlpha;                                                          // The alpha parameter for the simplex move on branch lengths
    private:
        bool                fixedTree;                                                          // Whether the tree topology is fixed
        double              currentPrior;                                                       // The current log prior probability
        double              lambda;                                                             // The rate parameter for the exponential prior on branch lengths
        double              oldPrior;                                                           // The old log prior probability
        int                 branchAcceptCount;                                                  // The number of accepted branch scale moves
        int                 branchCount;                                                        // The number of proposed branch scale moves
        TreeMoves           moveChoice;                                                         // The last move proposed
        int                 treeAcceptCount;                                                    // The number of accepted simplex moves on the whole tree
        int                 treeCount;                                                          // The number of proposed simplex moves on the whole tree
        TreeObject*         trees[2];                                                           // The active and inactive tree objects
};

#endif