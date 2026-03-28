#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <memory>
#include <string>

/**
 * @brief Enum for tracking the ID of the last proposed move
 */
enum TreeMoves{
    NO_MOVE = 1,
    BRANCH_PROPORTION_MOVE = 2,
    TREE_LENGTH_MOVE = 3
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
                            TreeParameter(Alignment& aln, std::string& newick, double p[2]);    // Construct a new tree from an alignment, a newick string, and the parameters for the prior
                            TreeParameter(TreeObject& tree, double p[2]);                       // Construct a new tree from a pre-made tree object and a rate parameter for the exponential prior
                            TreeParameter(const TreeParameter& t);                              // Deep copy constructor
                            ~TreeParameter();                                                   // Destructor

        double              getBranchAlpha() const {return tuningState->branchAlpha;}           // Get the shared alpha parameter for branch proposals
        double              getBranchRate() const;                                              // Get the acceptance rate for branch proposals
        double              getTreeDelta() const {return tuningState->treeDelta;}               // Get the shared delta parameter for tree-length proposals
        double              getTreeRate() const;                                                // Get the acceptance rate for proposals on the whole tree lengths
        double              lnPrior() override;                                                 // Return the log prior probability of the tree
        double              update();                                                           // Propose a random MCMC move
        std::string         writeNewick() const {return trees[0]->getNewick();}                 // Return the newick string of the active tree object
        TreeObject*         getTree() const {return trees[0];}                                  // Get the active tree object
        void                accept() override;                                                  // Accept changes proposed to the tree
        void                reject() override;                                                  // Reject changes proposed to the tree
        void                setBranchAlpha(double alpha) {tuningState->branchAlpha = alpha;}    // Set the shared alpha parameter for branch proposals
        void                setCountTuningEvents(bool shouldCount) {countTuningEvents = shouldCount;}
        void                setTreeDelta(double delta) {tuningState->treeDelta = delta;}        // Set the shared delta parameter for tree-length proposals
        void                shareTuningWith(TreeParameter& t) {tuningState = t.tuningState;}    // Share burn-in tuning state with another tree parameter
        void                tune() override;                                                    // Tune proposals that update the tree
    private:
        struct ProposalTuningStats {
            int acceptCount = 0;
            int count = 0;
        };
        
        struct TreeTuningState {
            double treeDelta = 1.0;
            double branchAlpha = 1000.0;
            ProposalTuningStats branchStats;
            ProposalTuningStats treeStats;
        };

        double              currentPrior = 0.0;                                                 // The current log prior probability
        double              priorParams[2];                                                     // The rate parameter for the exponential prior on branch lengths
        double              oldPrior = 0.0;                                                     // The old log prior probability
        int                 branchAcceptCount = 0;                                              // The number of accepted branch moves
        int                 branchCount = 0;                                                    // The number of proposed branch moves
        bool                countTuningEvents = true;                                           // Whether accepted/rejected moves should count toward shared tuning
        TreeMoves           moveChoice = TreeMoves::NO_MOVE;                                    // The last move proposed
        int                 treeAcceptCount = 0;                                                // The number of accepted tree length moves
        int                 treeCount = 0;                                                      // The number of proposed tree length moves
        std::shared_ptr<TreeTuningState> tuningState = std::make_shared<TreeTuningState>();    // Shared tuning parameters and posterior-chain burn-in counts
        TreeObject*         trees[2];                                                           // The active and inactive tree objects
};

#endif
