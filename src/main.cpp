#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/CodonMultiMatrix.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/analysis/Mcmc.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <algorithm>
#include <chrono>

int main(int argc, char* argv[]) {

    Settings settings(argc, argv);

    if(!settings.simulate){
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        RandomVariable& rng = RandomVariable::randomVariableInstance();
        Alignment aln(settings.nexusInput);
        std::cout << "Initializing model..." << std::endl;

        TreeParameter treeParam(&aln, settings.fixedTree, settings.treeLengthLambda);

        DirichletProcessPrior dpp(aln.getNumChar(), settings);

        CodonMultiMatrix rateMatrix(settings);

        Model model(&aln, &treeParam, &rateMatrix, &dpp);

        Mcmc myMCMC(&model, &treeParam, &rateMatrix, &dpp, settings);

        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
    }
    else {
        // TODO: Fix this when everything is a little more settled
        /*
        RandomVariable& rng = RandomVariable::randomVariableInstance();

        std::cout << "Generating Tree..." << std::endl;
        TreeObject tree(settings.numTaxa);
        std::vector<Node*> nodes = tree.getPostOrderSeq();
        for(Node* n : nodes) {
            if(n != tree.getRoot())
                tree.setBranchLength(n, Probability::Exponential::rv(&rng, settings.treeLengthLambda));
        }


        std::cout << "Tree: " << tree.getNewick() << std::endl;
        Node* root = tree.getRoot();

        std::cout << "Generating Root Sequence and Allocating Alignment Space..." << std::endl;

        int numChar = settings.numChar;

        int alignmentWidth = numChar * tree.getNumNodes();
        int* alignment = new int[alignmentWidth];
        std::fill(alignment, alignment + alignmentWidth, 0);

        std::vector<double> stationaryDist;
        std::vector<double> dist;
        for(int i = 0; i < 122; i++){
            stationaryDist.push_back(1);
            dist.push_back(1);
        }

        Probability::Dirichlet::rv(&rng, dist, stationaryDist);

        int* pR = alignment + (numChar * root->getIndex());
        for(int c = 0; c < numChar; c++){
            double randomChar = rng.uniformRv();
            double total = 0.0;
            for(int i = 0; i < 122; i++){
                total += stationaryDist[i];
                if(randomChar < total){
                    *pR = i;
                    pR += 1;
                    break;
                }
            }
        }

        std::vector<int> assignments = settings.assignmentVector;
        std::vector<double> omega1 = settings.omega1Vector;
        std::vector<double> omega2 = settings.omega2Vector;
        int numCats = omega1.size();

        if(assignments.size() == 0){
            std::cout << "No Assignments Provided! Generating a Random Assignment for " << numCats << " Categories." << std::endl;
            double increment = 1.0/numCats;;
            for(int c = 0; c < numChar; c++){
                double randomCat = rng.uniformRv();

                double total = 0.0;
                for(int i = 0; i < numCats; i++){
                    total += increment;
                    if(randomCat < total){
                        assignments.push_back(i);
                        std::cout << i << " ";
                        break;
                    }
                }
            }
            std::cout << std::endl;
        }

        std::cout << "Generating Rate Matrix and Transition Probabilities..." << std::endl;
        CodonMultiMatrix rateMatrix(settings, stationaryDist);
        TransitionProbability transProb(tree.getNumNodes(), numCats);

        std::vector<Node*>&  preOrderSeq = tree.getPostOrderSeq();
        std::reverse(preOrderSeq.begin(), preOrderSeq.end());

        for(int i = 0; i < numCats; i++){
            transProb.updateQ(rateMatrix.Q(omega1[i], omega2[i]), i);
        }

        for(Node* n : preOrderSeq){
            if(n != root){
                int nIndex = n->getIndex();
                double v = tree.getBranchLength(n);
                for(int i = 0; i < numCats; i++)
                    transProb.setProbs(0, i, nIndex, v);
            }
        }

        std::cout << "Simulated Pre-Order Traversal:" << std::endl;

        std::vector<const char*> codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

        for(Node* n : preOrderSeq){
            if(!n->getIsTip()){
                int nIndex = n->getIndex();
                int* pN = alignment + (numChar * n->getIndex());

                std::set<Node*>& nNeighbors = n->getNeighbors();
                for(Node* d : nNeighbors){
                    if(d != n->getAncestor()){
                        int dIndex = d->getIndex();
                        int* pD = alignment + (dIndex * numChar);
                        int* pNN = pN;

                        std::vector<Matrix<double>> transProbMatrices;
                        for(int cat = 0; cat < numCats; cat++)
                            transProbMatrices.push_back(*(transProb)(0, cat, dIndex));
                        for(int c = 0; c < numChar; c++){
                            Matrix<double> P = transProbMatrices[assignments[c]];

                            double randomChoice = rng.uniformRv();

                            double total = 0.0;
                            for(int i = 0; i < 122; i++){
                                total += P(*pNN, i);

                                if(randomChoice < total){
                                    *pD = i;
                                    pD += 1;
                                    pNN += 1;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            else {
                int nIndex = n->getIndex();
                int* pN = alignment + (numChar * n->getIndex());
                
                std::cout << ">" << n->getName() << std::endl;;
                for (int c = 0; c < numChar; c++){
                    int cIndex = *pN;
                    if(cIndex < 61)
                        std::cout << codons[cIndex];
                    else
                        std::cout << codons[cIndex - 61];
                    pN += 1;
                }
                std::cout << std::endl;
            }
        }

        delete [] alignment;
    */
    }
}