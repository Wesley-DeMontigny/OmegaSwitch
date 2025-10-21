#ifndef TREE_OBJECT_HPP
#define TREE_OBJECT_HPP
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

class Node;
class Alignment;
class RandomVariable;

/**
 * @brief 
 * 
 */
class TreeObject {

    public:
                                    TreeObject(void) = delete;
                                    TreeObject(int nt, bool rooted);
                                    TreeObject(Alignment* aln, bool rooted);
                                    TreeObject(const TreeObject& t);
                                    TreeObject(std::string newick, std::vector<std::string> taxaNames);
                                    ~TreeObject(void);
        TreeObject&                 operator=(const TreeObject& rhs);

        double                      getBranchLength(Node* n) const;
        int                         getNumNodes(){return nodes.size();}
        int                         getNumTaxa(){return numTaxa;}
        Node*                       getRoot() {return root;}
        std::map<Node*, double>     getBranchLengthMapping();
        std::string                 getNewick() const;
        std::vector<double>         getBranchLengths();
        std::vector<Node*>          getPostOrderSeq() {return postOrderSeq;}
        std::vector<Node*>          getTips();
        void                        initPostOrder(void);
        void                        passDown(Node* p, std::vector<Node*>& vec);
        void                        print(std::string header) const;
        void                        print(void) const;
        void                        setBranchLength(Node* n, double length);
        void                        updateAll();
        
    private:
        int                         getTaxonIndex(std::string token, std::vector<std::string> taxaNames);
        Node*                       addNode(void);
        std::vector<std::string>    parseNewickString(std::string newick);
        void                        clone(const TreeObject& t);
        void                        deleteAllNodes();
        void                        showNode(Node* p, int indent) const;
        void                        writeNode(Node* p, std::stringstream& strm) const;

        int                         numTaxa;
        Node*                       root;
        std::map<Node*, double>     branchLengths;
        std::vector<Node*>          nodes;
        std::vector<Node*>          postOrderSeq;
};

#endif
