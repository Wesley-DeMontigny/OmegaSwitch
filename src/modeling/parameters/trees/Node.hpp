#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include <string>
#include <unordered_map>

/**
 * @brief The nodes a phylogenetic tree consists of
 */
class Node{
    public:
                            Node();

        bool                getIsTip() const {return isTip;}
        bool                getNeedsCLUpdate(){return needsCLUpdate;}
        bool                getNeedsTPUpdate(){return needsTPUpdate;}
        int                 getIndex() const {return index;}
        int                 getOffset() const {return offset;}
        Node*               getAncestor() const {return ancestor;}
        static Node*        chooseNodeFromSet(std::set<Node*>& s);
        std::set<Node*>&    getNeighborRef() {return neighbors;}
        std::string         getName() const {return name;}
        void                addNeighbor(Node* n) {neighbors.insert(n);}
        void                removeAllNeighbors() {neighbors.clear();}
        void                removeNeighbor(Node* n) {neighbors.erase(n);}
        void                setAncestor(Node* a) {ancestor = a;}
        void                setIndex(int i) {index = i;}
        void                setIsTip(bool t) {isTip = t;}
        void                setName(std::string s) {name = s;}
        void                setNeedsCLUpdate(bool nU){needsCLUpdate = nU;}
        void                setNeedsTPUpdate(bool nU){needsTPUpdate = nU;}
        void                setOffset(int o){offset = o;}
    private:
        bool                isTip;
        bool                needsCLUpdate;
        bool                needsTPUpdate;
        int                 index;
        int                 offset;
        Node*               ancestor;
        std::set<Node*>     neighbors;
        std::string         name;
};


namespace std {
    /**
     * @brief Create a hash object so we can use unordered_maps with Node* keys
     */
    template <>
    struct hash<Node*> {
        std::size_t operator()(Node* k) const {
            return std::hash<int>()(k->getIndex());
        }
    };
}

#endif