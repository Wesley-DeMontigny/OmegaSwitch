#include "Node.hpp"
#include "core/RandomVariable.hpp"

Node::Node() : index(0), ancestor(nullptr), name(""), isTip(false), offset(0), needsCLUpdate(false), needsTPUpdate(false) {}

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
