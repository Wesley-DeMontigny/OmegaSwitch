#ifndef TREE_OBJECT_HPP
#define TREE_OBJECT_HPP
#include "Node.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <set>


class Alignment;
class RandomVariable;

/**
 * @brief The class containing the phylogenetic tree datastructure
 */
class TreeObject {
    public:
                                            TreeObject(void) = delete;
                                            TreeObject(int nt, bool rooted);                                        //
                                            TreeObject(Alignment& aln, bool rooted);                                //
                                            TreeObject(const TreeObject& t);                                        //
                                            TreeObject(std::string newick, std::vector<std::string> taxaNames);     //
                                            ~TreeObject();                                                          //
        TreeObject&                         operator=(const TreeObject& rhs);                                       //

        double                              getBranchLength(Node* n) const;                                         //
        double                              getTreeLength() const {return treeLength;}                              //
        int                                 getNumNodes() const {return nodes.size();}                              //
        int                                 getNumTaxa() const {return numTaxa;}                                    //
        Node*                               getRoot() const {return root;}                                          //
        std::unordered_map<Node*, double>   getBranchPropMapping() const;                                           //
        std::string                         getNewick() const;                                                      //
        std::vector<double>                 getBranchProportions() const;                                           //
        std::vector<Node*>                  getPostOrderSeq() const {return postOrderSeq;}                          //
        std::vector<Node*>                  getTips() const;                                                        //
        void                                initPostOrder();                                                        //
        void                                print(std::string header) const;                                        //
        void                                print() const;                                                          //
        void                                randomizeBranches();                                                    //
        void                                setBranchProportion(Node* n, double p);                                 //
        void                                setTreeLength(double length) {treeLength = length;};                    //
        void                                updateAll();                                                            //
        
    private:
        int                                 getTaxonIndex(std::string token, std::vector<std::string> taxaNames);   //
        Node*                               addNode(void);                                                          //
        std::vector<std::string>            parseNewickString(std::string newick);                                  //
        void                                clone(const TreeObject& t);                                             //
        void                                deleteAllNodes();                                                       //
        void                                passDown(Node* p, std::vector<Node*>& vec);                             //
        void                                showNode(Node* p, int indent) const;                                    //
        void                                writeNode(Node* p, std::stringstream& strm) const;                      //

        int                                 numTaxa;                                                                //
        Node*                               root;                                                                   //
        std::unordered_map<Node*, double>   branchProportions;                                                      //
        double                              treeLength = 6.0;                                                       //
        std::vector<Node*>                  nodes;                                                                  //
        std::vector<Node*>                  postOrderSeq;                                                           //
};

#endif
