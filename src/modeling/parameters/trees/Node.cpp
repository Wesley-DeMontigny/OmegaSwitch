#include "core/RandomVariable.hpp"
#include "Node.hpp"

/**
 * @brief Default constructor
 */
Node::Node() : index(0), ancestor(nullptr), name(""), isTip(false), offset(0), needsCLUpdate(false), needsTPUpdate(false) {}

/**
 * @brief Choose a node from a set of nodes. Useful for random sampling.
 */
Node* Node::chooseNodeFromSet(std::set<Node*>& s){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double rand = rng.uniformRv();
    int whichNode = (int)(s.size() * rand);
    int i = 0;
    for(Node* n : s){
        if(whichNode == i)
            return n;
        i++;
    }
    return nullptr;
}
