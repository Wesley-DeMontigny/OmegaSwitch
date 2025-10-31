#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "core/RandomVariable.hpp"
#include "Node.hpp"
#include "TreeObject.hpp"
#include <cmath>
#include <iostream>


TreeObject::TreeObject(int nt, bool rooted) : numTaxa(nt) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    root = addNode();
    root->setName("Root");
    int initTips = 0;
    if(rooted == false){
        initTips = 3;
    }
    else {
        initTips = 2;
    }

    for(int i = 0; i < initTips; i++) {
        Node* p = addNode();
        p->setName("Taxon_" + std::to_string(i+1));
        p->setIsTip(true);
        p->setIndex(i);
        p->addNeighbor(root);
        root->addNeighbor(p);
        p->setAncestor(root);
    }

    // build up the full tree by randomly adding branches to the existing tree
    for(int i = initTips; i < numTaxa; i++){
        Node* p = nullptr;
        do {
            p = nodes[(int)(rng.uniformRv() * nodes.size())];
        } while(p == root);

        Node* pAnc = p->getAncestor();
        Node* newTip = addNode();
        newTip->setIsTip(true);
        newTip->setIndex(i);
        newTip->setName("Taxon_" + std::to_string(i+1));
        Node* newAnc = addNode();

        p->removeNeighbor(pAnc);
        p->addNeighbor(newAnc);
        p->setAncestor(newAnc);

        newAnc->addNeighbor(p);
        newAnc->addNeighbor(newTip);
        newAnc->addNeighbor(pAnc);
        newAnc->setAncestor(pAnc);

        newTip->addNeighbor(newAnc);
        newTip->setAncestor(newAnc);

        pAnc->removeNeighbor(p);
        pAnc->addNeighbor(newAnc);
    }
        
    // initialize the down pass sequence
    initPostOrder();

    // Initialize branch lengths
    for (int i=0, n=(int)postOrderSeq.size(); i<n; i++) {
        Node* p = postOrderSeq[i];
        if (p->getAncestor() != nullptr)
                this->setBranchLength(p, rng.uniformRv());
    }

    // index the interior nodes (the tip nodes are indexed, above)
    int intIdx = numTaxa;
    for (int i=0, n=(int)postOrderSeq.size(); i<n; i++) {
        Node* p = postOrderSeq[i];
        if (p->getIsTip() == false)
            p->setIndex(intIdx++);
    }
}

//Generate a random tree and connect it to an alignment
TreeObject::TreeObject(Alignment& aln, bool rooted) : TreeObject(aln.getNumTaxa(), rooted) {

    std::vector<std::string> names = aln.getTaxaNames();
    for(Node* n : postOrderSeq){
        if(n->getIsTip())
            n->setName(names[n->getIndex()]);
    }
    
}

//Change to match index to taxon name
TreeObject::TreeObject(std::string newick, std::vector<std::string> taxaNames){
    std::vector<std::string> tokens = parseNewickString(newick);

    Node* p = nullptr;
    bool readingBranchLength = false;

    numTaxa = 0;

    for(std::string tok : tokens){
        if(tok == "("){
            Node* newNode = addNode();
            if(p == nullptr)
                root = newNode;
            else{
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                setBranchLength(newNode, 0.0);
            }

            p = newNode;
        }
        else if(tok == ")" || tok == ","){
            if(p->getAncestor() != nullptr)
                p = p->getAncestor();
            else
                Msg::error("Poorly formatted Newick! -P should not be null.");
        }
        else if(tok == ":"){
            readingBranchLength = true;
        }
        else if(tok == ";"){
            if(p != root)
                Msg::error("Poorly formatted Newick! Did not end at root.");
        }
        else{
            if(readingBranchLength){
                double x = atof(tok.c_str());
                setBranchLength(p, x);
            }
            else{
                //We need to trim the white space at the beginning and end of the token
                while(tok[0] == ' ')
                    tok.erase(0,1);
                while(tok[tok.size()-1] == ' ')
                    tok.erase(tok.size()-1);

                // There is a weird bug causing spaces elsewhere to not be parsed correctly? I am putting this here to fix that
                if(tok != ""){
                    Node* newNode = addNode();
                    p->addNeighbor(newNode);
                    newNode->addNeighbor(p);
                    newNode->setAncestor(p);
                    newNode->setName(tok);
                    newNode->setIsTip(true);
                    setBranchLength(newNode, 0.0);

                    int taxonIndex = getTaxonIndex(tok, taxaNames);
                    if(taxonIndex == -1)
                        Msg::error("Token '" + tok + "' is not in taxa names");
                    newNode->setIndex(taxonIndex);

                    p = newNode;
                    numTaxa++;
                }
            }
            readingBranchLength = false;
        }
    }
    initPostOrder();
    
    if(numTaxa != taxaNames.size())
        Msg::error("Taxa names do not match the size of the newick string.");

    int idx = numTaxa;
    for(Node* n : postOrderSeq){
        if(n->getIsTip() == false){
            n->setIndex(idx++);
        }
    }
}

TreeObject::TreeObject(const TreeObject& t){
    clone(t);
}

TreeObject::~TreeObject(void) {
    deleteAllNodes();
}

//Deep Copy Operation
TreeObject& TreeObject::operator=(const TreeObject& rhs){
    if(this == &rhs)
        return *this;
    
    clone(rhs);
    return *this;
}

Node* TreeObject::addNode(void) {

    Node* newNode = new Node;
    newNode->setOffset((int)nodes.size());
    nodes.push_back(newNode);
    return newNode;
}

void TreeObject::clone(const TreeObject& t){
    if(nodes.size() != t.nodes.size()){
        deleteAllNodes();
        for(int i = 0; i < t.nodes.size(); i++)
            addNode();
    }
    this->branchLengths.clear();

    this->numTaxa = t.numTaxa;
    this->root = this->nodes[t.root->getOffset()];

    for(int i = 0; i < t.nodes.size(); i++){
        Node* p = this->nodes[i];
        Node* q = t.nodes[i];
        p->setIndex(q->getIndex());
        p->setIsTip(q->getIsTip());
        p->setName(q->getName());
        p->setNeedsCLUpdate(q->getNeedsCLUpdate());
        p->setNeedsTPUpdate(q->getNeedsTPUpdate());

        p->removeAllNeighbors();
        std::set<Node*>& qNeighbors = q->getNeighborRef();
        for(Node* n : qNeighbors)
            p->addNeighbor(this->nodes[n->getOffset()]);

        if(q->getAncestor() != nullptr){
            Node* ancestor = this->nodes[q->getAncestor()->getOffset()];
            double bl = t.getBranchLength(q);
            p->setAncestor(ancestor);
            this->setBranchLength(p, bl);
        }
        else
            p->setAncestor(nullptr);
    }

    this->postOrderSeq.clear();
    for(int i = 0; i < t.postOrderSeq.size(); i++){
        Node* p = t.postOrderSeq[i];
        this->postOrderSeq.push_back(this->nodes[p->getOffset()]);
    }
}

void TreeObject::deleteAllNodes(){
    for (int i = 0; i < nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

void TreeObject::setBranchLength(Node* n, double length){
    auto it = branchLengths.find(n);

    if(it == branchLengths.end())
        branchLengths.insert(std::make_pair(n, length));
    else
        it->second = length;
}

double TreeObject::getBranchLength(Node* n) const{
    auto it = branchLengths.find(n);

    if(it == branchLengths.end())
        Msg::error("Couldn't find branch length of pair");
    return it->second;
}

std::vector<double> TreeObject::getBranchLengths() const {
    std::vector<double> returnVec;
    returnVec.reserve(branchLengths.size());

    for (auto s : branchLengths)
        returnVec.push_back(s.second);

    return returnVec;
}

std::unordered_map<Node*, double> TreeObject::getBranchLengthMapping() const {
    return branchLengths;
}

std::string TreeObject::getNewick() const{
    std::stringstream strm;
    writeNode(root, strm);
    return strm.str();
}

std::vector<Node*> TreeObject::getTips() const {
    std::vector<Node*> out;
    out.reserve(numTaxa);

    for(Node* n : nodes){
        if(n->getIsTip())
            out.push_back(n);
    }

    return out;
}

int TreeObject::getTaxonIndex(std::string token, std::vector<std::string> taxaNames){

    for(int i = 0, n = taxaNames.size(); i < n; i++){
        if(taxaNames[i] == token)
            return i;
    }

    return -1;
}

void TreeObject::initPostOrder(void) {
    postOrderSeq.clear();
    passDown(root, postOrderSeq);
}

std::vector<std::string> TreeObject::parseNewickString(std::string newick){
    std::vector<std::string> tokens;
    std::string str = "";
    for(int i = 0; i < newick.length(); i++){
        char c = newick[i];
        if(c == '(' || c == ')' || c == ',' || c == ':' || c == ';'){
            if(str != ""){
                tokens.push_back(str);
                str = "";
            }

            tokens.push_back(std::string(1, c));
        }
        else {
            str += std::string(1, c);
        }
    }

    return tokens;
}

void TreeObject::passDown(Node* p, std::vector<Node*>& vec) {

    if(p == nullptr)
        return;
    
    std::set<Node*>& pNeighbors = p->getNeighborRef();

    for(Node* n : pNeighbors){
        if(n != p->getAncestor())
            passDown(n, vec);
    }

    vec.push_back(p);
}

void TreeObject::print(std::string header) const{
    std::cout << header << std::endl;
    print();
}

void TreeObject::print(void) const{

    showNode(root, 0);
}

//Output nodes of the tree in a whitespace-indented format
void TreeObject::showNode(Node* p, int indent) const{

    if(p == nullptr)
        return;

    for(int i = 0; i < indent; i++)
        std::cout << " ";

    std::cout << p->getIndex() << " ( ";
    std::set<Node*>& pNeighbors = p->getNeighborRef();
    for (Node* d : pNeighbors)
        {
        if (d == p->getAncestor())
            std::cout << "a_";
        std::cout << d->getIndex() << " ";
        }
    std::cout << ") ";

    if(p->getAncestor() != nullptr)
        std::cout << this->getBranchLength(p) << " ";

    std::cout << p->getName();

    if (p == root)
        std::cout << " <- Root";
    std::cout << std::endl;

    for(Node* n : pNeighbors)
        {
        if(n != p->getAncestor())
            showNode(n, indent+3);
        }
}

void TreeObject::updateAll(){
    for(Node* n : nodes){
        if(n->getIsTip() == false)
            n->setNeedsCLUpdate(true);
        n->setNeedsTPUpdate(true);
    }
}

//For outputting a newick string
void TreeObject::writeNode(Node* p, std::stringstream& strm) const{
    if(p == nullptr)
        return;
    
    if(!p->getIsTip())
        strm << "(";
    else
        strm << p->getName() << "[&index=" << p->getIndex() << "]";

    std::set<Node*>& pDesc = p->getNeighborRef();
    bool foundFirst = false;
    for(Node* n : pDesc){
        if(n != p->getAncestor()){
            if(foundFirst)
                strm << ",";
            foundFirst = true;
            writeNode(n, strm);
        }
    }

    if(!p->getIsTip())
        strm << ")[&index=" << p->getIndex() << "]";
    if(p->getAncestor() != nullptr)
        strm << ":" << this->getBranchLength(p);
    else
        strm << ":0.0;";
}
