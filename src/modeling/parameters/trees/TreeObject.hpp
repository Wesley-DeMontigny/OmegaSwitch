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
                                            ~TreeObject(void);                                                      //
        TreeObject&                         operator=(const TreeObject& rhs);                                       //

        double                              getBranchLength(Node* n) const;                                         //
        int                                 getNumNodes() const {return nodes.size();}                              //
        int                                 getNumTaxa() const {return numTaxa;}                                    //
        Node*                               getRoot() const {return root;}                                          //
        std::unordered_map<Node*, double>   getBranchLengthMapping() const;                                         //
        std::string                         getNewick() const;                                                      //
        std::vector<double>                 getBranchLengths() const;                                               //
        std::vector<Node*>                  getPostOrderSeq() const {return postOrderSeq;}                          //
        std::vector<Node*>                  getTips() const;                                                        //
        void                                initPostOrder(void);                                                    //
        void                                print(std::string header) const;                                        //
        void                                print(void) const;                                                      //
        void                                setBranchLength(Node* n, double length);                                //
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
        std::unordered_map<Node*, double>   branchLengths;                                                          //
        std::vector<Node*>                  nodes;                                                                  //
        std::vector<Node*>                  postOrderSeq;                                                           //
};

#endif
