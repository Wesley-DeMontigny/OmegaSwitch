#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include <string>

class Node{
    public:
                            Node();
        void                addNeighbor(Node* n) {neighbors.insert(n);}
        static Node*        chooseNodeFromSet(std::set<Node*>& s);
        Node*               getAncestor() {return ancestor;}
        void                setAncestor(Node* a) {ancestor = a;}
        int                 getIndex() const {return index;}
        bool                getIsTip() const {return isTip;}
        std::string         getName() const {return name;}
        bool                getNeedsCLUpdate(){return needsCLUpdate;}
        bool                getNeedsTPUpdate(){return needsTPUpdate;}
        std::set<Node*>&    getNeighbors() {return neighbors;}
        int                 getOffset() const {return offset;} 
        void                removeNeighbor(Node* n) {neighbors.erase(n);}
        void                removeAllNeighbors() {neighbors.clear();}
        void                setIndex(int i) {index = i;}
        void                setIsTip(bool t) {isTip = t;}
        void                setName(std::string s) {name = s;}
        void                setNeedsCLUpdate(bool nU){needsCLUpdate = nU;}
        void                setNeedsTPUpdate(bool nU){needsTPUpdate = nU;}
        void                setOffset(int o){offset = o;}
    private:
        Node*               ancestor;
        int                 index;
        bool                isTip;
        std::string         name;
        bool                needsCLUpdate;
        bool                needsTPUpdate;
        std::set<Node*>     neighbors;
        int                 offset;
};

#endif