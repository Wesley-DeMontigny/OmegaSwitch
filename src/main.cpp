#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "core/Matrix.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/parameters/DPPMatrix.hpp"
#include "modeling/model/DPPModel.hpp"
#include "modeling/analysis/DPPMcmc.hpp"
#include "modeling/parameters/M0Matrix.hpp"
#include "modeling/model/M0Model.hpp"
#include "modeling/analysis/M0Mcmc.hpp"
#include "modeling/parameters/M3S2Matrix.hpp"
#include "modeling/model/M3S2Model.hpp"
#include "modeling/analysis/M3S2Mcmc.hpp"
#include "modeling/parameters/SBMatrix.hpp"
#include "modeling/model/SBModel.hpp"
#include "modeling/analysis/SBMcmc.hpp"
#include "modeling/parameters/RJMatrix.hpp"
#include "modeling/model/RJModel.hpp"
#include "modeling/analysis/RJMcmc.hpp"
#include "modeling/parameters/RJDPPMatrix.hpp"
#include "modeling/model/RJDPPModel.hpp"
#include "modeling/analysis/RJDPPMcmc.hpp"
#include "modeling/parameters/RJDirichletProcessPrior.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <algorithm>
#include <chrono>

#ifdef __AVX2__
#pragma message("Optimizing using AVX2")
#elif defined(__ARM_NEON__)
#pragma message("Optimizing using ARM NEON")
#else
#pragma message("No CPU optimizations available")
#endif

void inference(Settings& settings, Alignment& aln, TreeParameter& treeParam, bool disableBayesOpt){
    if(settings.M0){
        std::cout << "Initializing the M0 model..." << std::endl;

        M0Matrix rateMatrix(settings);

        M0Model model(settings, &aln, &treeParam, &rateMatrix);

        M0Mcmc myMCMC(&model, &treeParam, &rateMatrix, settings, disableBayesOpt);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else if(settings.M3S2){
        std::cout << "Initializing the M3S2 model..." << std::endl;

        M3S2Matrix rateMatrix(settings);

        M3S2Model model(settings, &aln, &treeParam, &rateMatrix);

        M3S2Mcmc myMCMC(&model, &treeParam, &rateMatrix, settings, disableBayesOpt);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else if(settings.SB){
        std::cout << "Initializing the SB model..." << std::endl;

        SBMatrix rateMatrix(settings);

        SBModel model(settings, &aln, &treeParam, &rateMatrix);

        SBMcmc myMCMC(&model, &treeParam, &rateMatrix, settings, disableBayesOpt);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else if(settings.RJDPP){
        std::cout << "Initializing the Reversible Jump DPP model..." << std::endl;

        RJDPPMatrix rateMatrix(settings);

        RJDirichletProcessPrior dpp(&rateMatrix, aln.getNumChar(), settings);

        RJDPPModel model(settings, &aln, &treeParam, &rateMatrix, &dpp);

        RJDPPMcmc myMCMC(&model, &treeParam, &rateMatrix, &dpp, settings, disableBayesOpt);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else if(settings.RJ){
        std::cout << "Initializing the Reversible Jump model..." << std::endl;

        RJMatrix rateMatrix(settings);

        RJModel model(settings, &aln, &treeParam, &rateMatrix);

        RJMcmc myMCMC(&model, &treeParam, &rateMatrix, settings, disableBayesOpt);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else {
        std::cout << "Initializing the DPP model..." << std::endl;
        
        DirichletProcessPrior dpp(aln.getNumChar(), settings);

        DPPMatrix rateMatrix(settings);

        DPPModel model(settings, &aln, &treeParam, &rateMatrix, &dpp);

        DPPMcmc myMCMC(&model, &treeParam, &rateMatrix, &dpp, settings, disableBayesOpt);

        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
}

int randomMultinomial(RandomVariable& rng, const std::vector<double>& weights){
    int iter = 0;
    do{
        double choice = rng.uniformRv();

        double cumulative = 0.0;
        for(int i = 0; i < weights.size(); i++){
            cumulative += weights[i];
            if(choice < cumulative){
                return i;
            }
        }
        iter++;
    }
    while(iter < 100);
    Msg::error("Failed 100 Multinomial draws!");
    return -1;
}

int randomTransition(RandomVariable& rng, const int ancestralState, const Matrix<double>& transitionProbs){
    int iter = 0;
    do{
        double choice = rng.uniformRv();

        double cumulative = 0.0;
        for(int i = 0; i < transitionProbs.dim1(); i++){
            cumulative += transitionProbs(ancestralState, i);
            if(choice < cumulative){
                return i;
            }
        }
        iter++;
    }
    while(iter < 100);

    Msg::error("Failed 100 transition draws!");
    return -1;
}

int main(int argc, char* argv[]) {

    Settings settings(argc, argv);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(settings.simulateDPP == false && settings.simulateM0 == false && settings.simulateM3S2 == false && settings.simulateSB == false && settings.simulateRJ == false && settings.simulateRJDPP == false){
        Alignment aln(settings.nexusInput);
        TreeParameter treeParam(&aln, settings.fixedTree, settings.treeLengthLambda);
        inference(settings, aln, treeParam, false);
    }
    else{
        int numSites = 100;
        int numTaxa = 25;
        std::string mcmcOutput = settings.mcmcOutput;
        std::string tipsOutput = settings.tipsOutput;
        std::string dppOutput = settings.dppOutput;
        std::string ancestralStatesOutput = settings.ancestralStatesOutput;
        std::string treeOutput = settings.treeOutput;
        std::string branchOutput = settings.branchOutput;
        
        std::string modelClass = "";
        if(settings.simulateDPP)
            modelClass = "DPP";
        else if(settings.simulateM0)
            modelClass = "M0";
        else if(settings.simulateM3S2)
            modelClass = "M3S2";
        else if(settings.simulateSB)
            modelClass = "SB";
        else if(settings.simulateRJ)
            modelClass = "RJ";
        else if(settings.simulateRJDPP)
            modelClass = "RJDPP";

        std::cout << "Doing a simulation analysis of " << modelClass << " with " << settings.numSimulations << " simulations of " << numTaxa << " taxa trees with " << numSites << " sites." << std::endl;
        for(int i = 0; i < settings.numSimulations; i++){
            std::cout << "Starting simulation " << i << "..." << std::endl;

            TreeObject tree(numTaxa, true);
            std::vector<Node*> preOrderSeq = tree.getPostOrderSeq();
            std::reverse(preOrderSeq.begin(), preOrderSeq.end());
            Node* root = tree.getRoot();
            int numNodes = tree.getNumNodes();

            for(Node* n : preOrderSeq) {
                if(n != root) {
                    tree.setBranchLength(n, Probability::Exponential::rv(&rng, settings.treeLengthLambda));
                }
            }

            int* sites = new int[numSites * numNodes];
            for(int c = 0; c < numSites * numNodes; c++){
                sites[c] = 0;
            }
            double* truedNdS = new double[numTaxa * numSites];
            for(int t = 0; t < numTaxa * numSites; t++){
                truedNdS[t] = 0.0;
            }
            int* tipSites = new int[numTaxa * numSites];
            for(int t = 0; t < numTaxa * numSites; t++){
                tipSites[t] = 0;
            }

            if(settings.simulateDPP){
                TransitionProbability transProb(numNodes, 122);
                DPPMatrix rateMatrix(settings);
                DirichletProcessPrior dpp(numSites, settings);
                std::vector<double> stationary = rateMatrix.getStationary();
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                std::vector<Category> categories = dpp.getCategories();
                std::vector<int> assignments = dpp.getAssignments();
                transProb.allocateQ(categories.size());
                std::vector<double> dNdS1;
                std::vector<double> dNdS2;
                for(int cat = 0; cat < categories.size(); cat++){
                    transProb.updateQ(rateMatrix.Q(categories[cat].omega1, categories[cat].omega2), cat);
                    auto dNdS = rateMatrix.dNdS(categories[cat].omega1, categories[cat].omega2);
                    dNdS1.push_back(std::get<0>(dNdS));
                    dNdS2.push_back(std::get<1>(dNdS));
                }

                for(Node* n : preOrderSeq){
                    if(n != root){
                        for(int cat = 0; cat < categories.size(); cat++){
                            transProb.setProbs(0, cat, n->getIndex(), tree.getBranchLength(n));
                        }
                    }
                }

                for(int c = 0; c < numSites; c++){
                    int cat = assignments[c];
                    for(Node* n : preOrderSeq){
                        int nIndex = n->getIndex();
                        if(n == root){
                            sites[c*numNodes + nIndex] = randomMultinomial(rng, stationary);
                        }
                        int ancestralState = sites[c*numNodes + nIndex];
                        if(n->getIsTip() == false){
                            for(Node* d : n->getNeighbors()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    Matrix<double> P = transProb(0, cat, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            if(sites[c*numNodes + nIndex] < 61){
                                truedNdS[c*numTaxa + nIndex] = dNdS1[cat];
                            }
                            else {
                                truedNdS[c*numTaxa + nIndex] = dNdS2[cat];
                            }
                            tipSites[c*numTaxa + nIndex] = sites[c*numNodes + nIndex];
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "K\tR\tCategoryCount";
                for(int element = 0; element < 61; element++){
                    fs << "\tStationary[" << element << "]";
                }
                fs << "\tTree";
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\tTaxon_" << taxon + 1 << "[" << site << "]";
                    }
                }
                fs << std::endl;
                fs << rateMatrix.getK() << "\t" << rateMatrix.getR() << "\t" << categories.size();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << stationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs.close();
            }
            else if(settings.simulateRJDPP){
                TransitionProbability transProb(numNodes, 183);
                RJDPPMatrix rateMatrix(settings);
                int activeOmegas = (int)(rng.uniformRv() * 3.0) + 1;
                RJDirichletProcessPrior dpp(&rateMatrix, numSites, settings);
                rateMatrix.refreshQBackground(activeOmegas);
                std::vector<double> stationary = rateMatrix.getStationary(activeOmegas);
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                std::vector<RJCategory> categories = dpp.getCategories();
                std::vector<int> assignments = dpp.getAssignments();
                transProb.allocateQ(categories.size());
                std::vector<double> dNdS1;
                std::vector<double> dNdS2;
                std::vector<double> dNdS3;
                for(int cat = 0; cat < categories.size(); cat++){
                    transProb.updateQ(rateMatrix.Q(categories[cat].omega1, categories[cat].omega2, categories[cat].omega3), cat);
                    auto dNdS = rateMatrix.dNdS(categories[cat].omega1, categories[cat].omega2, categories[cat].omega3);
                    dNdS1.push_back(std::get<0>(dNdS));
                    dNdS2.push_back(std::get<1>(dNdS));
                    dNdS3.push_back(std::get<2>(dNdS));
                }

                for(Node* n : preOrderSeq){
                    if(n != root){
                        for(int cat = 0; cat < categories.size(); cat++){
                            transProb.setProbs(0, cat, n->getIndex(), activeOmegas*61, tree.getBranchLength(n));
                        }
                    }
                }

                for(int c = 0; c < numSites; c++){
                    int cat = assignments[c];
                    for(Node* n : preOrderSeq){
                        int nIndex = n->getIndex();
                        if(n == root){
                            sites[c*numNodes + nIndex] = randomMultinomial(rng, stationary);
                        }
                        int ancestralState = sites[c*numNodes + nIndex];
                        if(n->getIsTip() == false){
                            for(Node* d : n->getNeighbors()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    Matrix<double> P = transProb(0, cat, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            if(sites[c*numNodes + nIndex] < 61){
                                truedNdS[c*numTaxa + nIndex] = dNdS1[cat];
                            }
                            else if(sites[c*numNodes + nIndex] < 122){
                                truedNdS[c*numTaxa + nIndex] = dNdS2[cat];
                            }
                            else {
                                truedNdS[c*numTaxa + nIndex] = dNdS3[cat];
                            }
                            tipSites[c*numTaxa + nIndex] = sites[c*numNodes + nIndex];
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "OmegaCount\tK\tR\tCategoryCount";
                for(int element = 0; element < 61; element++){
                    fs << "\tStationary[" << element << "]";
                }
                fs << "\tTree";
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\tTaxon_" << taxon + 1 << "[" << site << "]";
                    }
                }
                fs << std::endl;
                fs << activeOmegas << "\t" << rateMatrix.getK() << "\t" << rateMatrix.getR() << "\t" << categories.size();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << stationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs.close();
            }
            else if(settings.simulateM0){
                TransitionProbability transProb(numNodes, 61);
                M0Matrix rateMatrix(settings);
                std::vector<double> stationary = rateMatrix.getStationary();
                transProb.updateQ(rateMatrix.Q(), 0);
                double dNdS = rateMatrix.dNdS();

                for(Node* n : preOrderSeq){
                    if(n != root){
                        transProb.setProbs(0, 0, n->getIndex(), tree.getBranchLength(n));
                    }
                }

                for(int c = 0; c < numSites; c++){
                    for(Node* n : preOrderSeq){
                        int nIndex = n->getIndex();
                        if(n == root){
                            sites[c*numNodes + nIndex] = randomMultinomial(rng, stationary);
                        }
                        int ancestralState = sites[c*numNodes + nIndex];
                        if(n->getIsTip() == false){
                            for(Node* d : n->getNeighbors()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    Matrix<double> P = transProb(0, 0, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            truedNdS[c*numTaxa + nIndex] = dNdS;
                            tipSites[c*numTaxa + nIndex] = sites[c*numNodes + nIndex];
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "K\tOmega";
                for(int element = 0; element < 61; element++){
                    fs << "\tStationary[" << element << "]";
                }
                fs << "\tTree";
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\tTaxon_" << taxon + 1 << "[" << site << "]";
                    }
                }
                fs << std::endl;
                fs << rateMatrix.getK() << "\t" << rateMatrix.getOmega();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << stationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs.close();
            }
            else if(settings.simulateM3S2){
                TransitionProbability transProb(numNodes, 183);
                M3S2Matrix rateMatrix(settings);
                std::vector<double> stationary = rateMatrix.getStationary();
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                transProb.updateQ(rateMatrix.Q(), 0);
                auto dNdS = rateMatrix.dNdS();

                for(Node* n : preOrderSeq){
                    if(n != root){
                        transProb.setProbs(0, 0, n->getIndex(), tree.getBranchLength(n));
                    }
                }

                for(int c = 0; c < numSites; c++){
                    for(Node* n : preOrderSeq){
                        int nIndex = n->getIndex();
                        if(n == root){
                            sites[c*numNodes + nIndex] = randomMultinomial(rng, stationary);
                        }
                        int ancestralState = sites[c*numNodes + nIndex];
                        if(n->getIsTip() == false){
                            for(Node* d : n->getNeighbors()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    Matrix<double> P = transProb(0, 0, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            if(sites[c*numNodes + nIndex] < 61){
                                truedNdS[c*numTaxa + nIndex] = std::get<0>(dNdS);
                            }
                            else if(sites[c*numNodes + nIndex] < 122){
                                truedNdS[c*numTaxa + nIndex] = std::get<1>(dNdS);
                            }
                            else {
                                truedNdS[c*numTaxa + nIndex] = std::get<2>(dNdS);
                            }
                            tipSites[c*numTaxa + nIndex] = sites[c*numNodes + nIndex];
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "K\tOmega\tOmegaIncrement1\tOmegaIncrement2\tGamma\tR[1]\tR[2]";
                for(int element = 0; element < 61; element++){
                    fs << "\tStationary[" << element << "]";
                }
                fs << "\tTree";
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\tTaxon_" << taxon + 1 << "[" << site << "]";
                    }
                }
                fs << std::endl;
                fs << rateMatrix.getK() << "\t" << rateMatrix.getOmega1() 
                   << "\t" << rateMatrix.getOmega2() << "\t" << rateMatrix.getOmega3()
                   << "\t" << rateMatrix.getGamma() << "\t" << rateMatrix.getR1()
                   << "\t" << rateMatrix.getR2();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << rawStationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs.close();
            }
            else if(settings.simulateRJ){
                TransitionProbability transProb(numNodes, 183);
                RJMatrix rateMatrix(settings);
                int activeOmegas = (int)(rng.uniformRv() * 3.0) + 1;
                rateMatrix.setActiveOmegas(activeOmegas);
                std::vector<double> stationary = rateMatrix.getStationary();
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                transProb.updateQ(rateMatrix.Q(), 0);
                auto dNdS = rateMatrix.dNdS();

                for(Node* n : preOrderSeq){
                    if(n != root){
                        transProb.setProbs(0, 0, n->getIndex(), activeOmegas*61, tree.getBranchLength(n));
                    }
                }

                for(int c = 0; c < numSites; c++){
                    for(Node* n : preOrderSeq){
                        int nIndex = n->getIndex();
                        if(n == root){
                            sites[c*numNodes + nIndex] = randomMultinomial(rng, stationary);
                        }
                        int ancestralState = sites[c*numNodes + nIndex];
                        if(n->getIsTip() == false){
                            for(Node* d : n->getNeighbors()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    Matrix<double> P = transProb(0, 0, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            if(sites[c*numNodes + nIndex] < 61){
                                truedNdS[c*numTaxa + nIndex] = std::get<0>(dNdS);
                            }
                            else if(sites[c*numNodes + nIndex] < 122){
                                truedNdS[c*numTaxa + nIndex] = std::get<1>(dNdS);
                            }
                            else {
                                truedNdS[c*numTaxa + nIndex] = std::get<2>(dNdS);
                            }
                            tipSites[c*numTaxa + nIndex] = sites[c*numNodes + nIndex];
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "OmegaCount\tK\tOmega\tOmegaIncrement1\tOmegaIncrement2\tR";
                for(int element = 0; element < 61; element++){
                    fs << "\tStationary[" << element << "]";
                }
                fs << "\tTree";
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\tTaxon_" << taxon + 1 << "[" << site << "]";
                    }
                }
                fs << std::endl;
                fs << activeOmegas << "\t" << rateMatrix.getK() << "\t" << rateMatrix.getOmega1() 
                   << "\t" << rateMatrix.getOmega2() << "\t" << rateMatrix.getOmega3()
                   << "\t" << rateMatrix.getR();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << rawStationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs.close();
            }
            else if(settings.simulateSB){
                int numStates = settings.truncation * 61;
                TransitionProbability transProb(numNodes, numStates);
                SBMatrix rateMatrix(settings);
                std::vector<double> stationary = rateMatrix.getStationary();
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                transProb.updateQ(rateMatrix.Q(), 0);
                auto dNdS = rateMatrix.dNdS();

                for(Node* n : preOrderSeq){
                    if(n != root){
                        transProb.setProbs(0, 0, n->getIndex(), tree.getBranchLength(n));
                    }
                }

                for(int c = 0; c < numSites; c++){
                    for(Node* n : preOrderSeq){
                        int nIndex = n->getIndex();
                        if(n == root){
                            sites[c*numNodes + nIndex] = randomMultinomial(rng, stationary);
                        }
                        int ancestralState = sites[c*numNodes + nIndex];
                        if(n->getIsTip() == false){
                            for(Node* d : n->getNeighbors()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    Matrix<double> P = transProb(0, 0, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            int dNdSIndex = sites[c*numNodes + nIndex] / 61;
                            truedNdS[c*numTaxa + nIndex] = dNdS[dNdSIndex];
                            tipSites[c*numTaxa + nIndex] = sites[c*numNodes + nIndex];
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "K\tR";
                for(int o = 0; o < settings.truncation; o++){
                    fs << "\tOmega[" << o << "]";
                    fs << "\tProportion[" << o << "]";
                }
                for(int element = 0; element < 61; element++){
                    fs << "\tStationary[" << element << "]";
                }
                fs << "\tTree";
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\tTaxon_" << taxon + 1 << "[" << site << "]";
                    }
                }
                fs << std::endl;
                fs << rateMatrix.getK() << "\t" << rateMatrix.getR();
                std::vector<double> omegas = rateMatrix.getOmegas();
                std::vector<double> proportions = rateMatrix.getProportions();
                for(int o = 0; o < settings.truncation; o++){
                    fs << "\t" << omegas[o];
                    fs << "\t" << proportions[o];
                }
                for(int element = 0; element < 61; element++){
                    fs << "\t" << rawStationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs.close();
            }

            std::cout << "Simulation " << i << " is complete. Starting inference..." << std::endl;
            Alignment aln(tipSites, numSites, numTaxa);
            delete [] sites;
            delete [] tipSites;
            delete [] truedNdS;

            //Change the name I am logging to according to the number simulation it is
            if(!settings.sequentialTuningSim){
                TreeParameter treeParam(tree, settings.treeLengthLambda);
                if(settings.treeOutput != "")
                    settings.treeOutput = treeOutput + std::to_string(i);
                if(settings.mcmcOutput != "")
                    settings.mcmcOutput = mcmcOutput + std::to_string(i);
                if(settings.dppOutput != "")
                    settings.dppOutput = dppOutput + std::to_string(i);
                if(settings.tipsOutput != "")
                    settings.tipsOutput = tipsOutput + std::to_string(i);
                if(settings.ancestralStatesOutput != "")
                    settings.ancestralStatesOutput = ancestralStatesOutput + std::to_string(i);
                if(settings.branchOutput != "")
                    settings.branchOutput = branchOutput + std::to_string(i);
                
                inference(settings, aln, treeParam, false);
            }
            else{
                if(settings.treeOutput != "")
                    settings.treeOutput = treeOutput + std::to_string(i) + "_Bayes";
                if(settings.mcmcOutput != "")
                    settings.mcmcOutput = mcmcOutput + std::to_string(i) + "_Bayes";
                if(settings.dppOutput != "")
                    settings.dppOutput = dppOutput + std::to_string(i) + "_Bayes";
                if(settings.tipsOutput != "")
                    settings.tipsOutput = tipsOutput + std::to_string(i) + "_Bayes";
                if(settings.ancestralStatesOutput != "")
                    settings.ancestralStatesOutput = ancestralStatesOutput + std::to_string(i) + "_Bayes";
                if(settings.branchOutput != "")
                    settings.branchOutput = branchOutput + std::to_string(i) + "_Bayes";
                
                TreeParameter treeParamA(tree, settings.treeLengthLambda);
                inference(settings, aln, treeParamA, false);

                if(settings.treeOutput != "")
                    settings.treeOutput = treeOutput + std::to_string(i) + "_Classic";
                if(settings.mcmcOutput != "")
                    settings.mcmcOutput = mcmcOutput + std::to_string(i) + "_Classic";
                if(settings.dppOutput != "")
                    settings.dppOutput = dppOutput + std::to_string(i) + "_Classic";
                if(settings.tipsOutput != "")
                    settings.tipsOutput = tipsOutput + std::to_string(i) + "_Classic";
                if(settings.ancestralStatesOutput != "")
                    settings.ancestralStatesOutput = ancestralStatesOutput + std::to_string(i) + "_Classic";
                if(settings.branchOutput != "")
                    settings.branchOutput = branchOutput + std::to_string(i) + "_Classic";
                
                TreeParameter treeParamB(tree, settings.treeLengthLambda);
                inference(settings, aln, treeParamB, true);
            }
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
}