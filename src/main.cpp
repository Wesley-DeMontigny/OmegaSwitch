#include "misc/Alignment.hpp"
#include "misc/Matrix.hpp"
#include "misc/Probability.hpp"
#include "misc/RandomVariable.hpp"
#include "misc/Settings.hpp"
#include "modeling/analysis/MCMC.hpp"
#include "modeling/analysis/Move.hpp"
#include "modeling/model/M0Model.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/model/DPCMMModel.hpp"
#include "modeling/model/CMMModel.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/parameters/rate_matrices/M0Matrix.hpp"
#include "modeling/parameters/rate_matrices/DPCMMMatrix.hpp"
#include "modeling/parameters/rate_matrices/CMMMatrix.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "ncl/nxscharactersblock.h"
#include <algorithm>
#include <chrono>
#include <taskflow/taskflow.hpp>

#ifdef __AVX2__
#pragma message("Optimizing using AVX2")
#elif defined(__ARM_NEON__)
#pragma message("Optimizing using ARM NEON")
#else
#pragma message("No CPU optimizations available")
#endif

void inference(Settings& settings, Alignment& aln, TreeParameter& treeParam, bool disableBayesOpt, tf::Executor& executor){
    bool useMCMCMC = settings.mcmcmcBeta < 1.0;
    double alphaLength = (settings.treeLengthMean * settings.treeLengthMean) / (settings.treeLengthSD * settings.treeLengthSD);
    double betaLength = settings.treeLengthMean / (settings.treeLengthSD * settings.treeLengthSD);
    double treeLengthParams[2] = {alphaLength, betaLength};

    if(settings.CMM){
        std::cout << "Initializing the CMM model..." << std::endl;

        CMMMatrix rateMatrix(settings);

        CMMModel model(&settings, &aln, &treeParam, &rateMatrix, executor);

        std::vector<Move> moves;
        moves.emplace_back(Move{
            settings.treeWeight,
            static_cast<int>(aln.getNumTaxa() / 2.0),
            [&treeParam]() {return treeParam.update();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.omegaWeight,
            5,
            [&rateMatrix]() {return rateMatrix.updateOmega();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.kWeight,
            3,
            [&rateMatrix]() {return rateMatrix.updateK();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.stationaryWeight,
            10,
            [&rateMatrix]() {return rateMatrix.updateStationary();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.rjWeight,
            10,
            [&rateMatrix]() {return rateMatrix.updateActiveOmegas();},
            [&rateMatrix]() {return !rateMatrix.isRegimeCountFixed();}
        });
        moves.emplace_back(Move{
            settings.rWeight,
            5,
            [&rateMatrix]() {return rateMatrix.updateR();},
            [&rateMatrix]() {return rateMatrix.getActiveOmegas() > 1;}
        });

        if(useMCMCMC){
            TreeParameter temperedTree(treeParam);
            temperedTree.shareTuningWith(treeParam);
            CMMMatrix temperedRateMatrix(settings);
            temperedRateMatrix.shareTuningWith(rateMatrix);
            CMMModel temperedModel(&settings, &aln, &temperedTree, &temperedRateMatrix, executor);

            std::vector<Move> temperedMoves;
            temperedMoves.emplace_back(Move{
                settings.treeWeight,
                static_cast<int>(aln.getNumTaxa() / 2.0),
                [&temperedTree]() {return temperedTree.update();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.omegaWeight,
                5,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateOmega();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.kWeight,
                3,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateK();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.stationaryWeight,
                10,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateStationary();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.rjWeight,
                10,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateActiveOmegas();},
                [&temperedRateMatrix]() {return !temperedRateMatrix.isRegimeCountFixed();}
            });
            temperedMoves.emplace_back(Move{
                settings.rWeight,
                5,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateR();},
                [&temperedRateMatrix]() {return temperedRateMatrix.getActiveOmegas() > 1;}
            });

            MCMC mcmc(&model, moves, &temperedModel, temperedMoves, settings, disableBayesOpt);
            std::cout << "Starting MCMCMC with beta " << settings.mcmcmcBeta
                      << " and swap frequency " << settings.mcmcmcSwapFrequency << "..." << std::endl;
            mcmc.burnin();
            mcmc.run();
        }
        else{
            MCMC mcmc(&model, moves, settings, disableBayesOpt);
            
            std::cout << "Starting MCMC..." << std::endl;
            mcmc.burnin();
            mcmc.run();
        }
    }
    else if(settings.DPCMM){
        std::cout << "Initializing the DPCMM model..." << std::endl;

        DPCMMMatrix rateMatrix(settings, aln.getNumChar());

        DPCMMModel model(&settings, &aln, &treeParam, &rateMatrix, executor);

        std::vector<Move> moves;
        moves.emplace_back(Move{
            settings.treeWeight,
            static_cast<int>(aln.getNumTaxa() / 2.0),
            [&treeParam]() {return treeParam.update();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.omegaWeight,
            10,
            [&rateMatrix]() {return rateMatrix.updateOmega();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.kWeight,
            5,
            [&rateMatrix]() {return rateMatrix.updateK();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.stationaryWeight,
            10,
            [&rateMatrix]() {return rateMatrix.updateStationary();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.rjWeight,
            10,
            [&rateMatrix]() {return rateMatrix.updateActiveOmegas();},
            [&rateMatrix]() {return !rateMatrix.isRegimeCountFixed();}
        });
        moves.emplace_back(Move{
            settings.rWeight,
            5,
            [&rateMatrix]() {return rateMatrix.updateR();},
            [&rateMatrix]() {return rateMatrix.getActiveOmegas() > 1;}
        });
        #ifndef SAMPLE_PRIOR
        moves.emplace_back(Move{
            settings.dppWeight,
            1,
            [&model]() {return model.updateDPP();},
            [&rateMatrix]() {return !rateMatrix.isAssignmentFixed();}
        });
        #endif

        if(useMCMCMC){
            TreeParameter temperedTree(treeParam);
            temperedTree.shareTuningWith(treeParam);
            DPCMMMatrix temperedRateMatrix(settings, aln.getNumChar());
            temperedRateMatrix.shareTuningWith(rateMatrix);
            DPCMMModel temperedModel(&settings, &aln, &temperedTree, &temperedRateMatrix, executor);

            std::vector<Move> temperedMoves;
            temperedMoves.emplace_back(Move{
                settings.treeWeight,
                static_cast<int>(aln.getNumTaxa() / 2.0),
                [&temperedTree]() {return temperedTree.update();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.omegaWeight,
                10,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateOmega();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.kWeight,
                5,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateK();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.stationaryWeight,
                10,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateStationary();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.rjWeight,
                10,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateActiveOmegas();},
                [&temperedRateMatrix]() {return !temperedRateMatrix.isRegimeCountFixed();}
            });
            temperedMoves.emplace_back(Move{
                settings.rWeight,
                5,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateR();},
                [&temperedRateMatrix]() {return temperedRateMatrix.getActiveOmegas() > 1;}
            });
            #ifndef SAMPLE_PRIOR
            temperedMoves.emplace_back(Move{
                settings.dppWeight,
                1,
                [&temperedModel]() {return temperedModel.updateDPP();},
                [&temperedRateMatrix]() {return !temperedRateMatrix.isAssignmentFixed();}
            });
            #endif

            MCMC mcmc(&model, moves, &temperedModel, temperedMoves, settings, disableBayesOpt);
            std::cout << "Starting MCMCMC with beta " << settings.mcmcmcBeta
                      << " and swap frequency " << settings.mcmcmcSwapFrequency << "..." << std::endl;
            mcmc.burnin();
            mcmc.run();
        }
        else{
            MCMC mcmc(&model, moves, settings, disableBayesOpt);
            
            std::cout << "Starting MCMC..." << std::endl;
            mcmc.burnin();
            mcmc.run();
        }
    }
    else { // Default to the M0 model
        std::cout << "Initializing the M0 model..." << std::endl;

        M0Matrix rateMatrix(settings);

        M0Model model(&settings, &aln, &treeParam, &rateMatrix, executor);

        std::vector<Move> moves;
        moves.emplace_back(Move{
            settings.treeWeight,
            aln.getNumTaxa(),
            [&treeParam]() {return treeParam.update();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.omegaWeight,
            5,
            [&rateMatrix]() {return rateMatrix.updateOmega();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.kWeight,
            3,
            [&rateMatrix]() {return rateMatrix.updateK();},
            []() {return true;}
        });
        moves.emplace_back(Move{
            settings.stationaryWeight,
            10,
            [&rateMatrix]() {return rateMatrix.updateStationary();},
            []() {return true;}
        });

        if(useMCMCMC){
            TreeParameter temperedTree(treeParam);
            temperedTree.shareTuningWith(treeParam);
            M0Matrix temperedRateMatrix(settings);
            temperedRateMatrix.shareTuningWith(rateMatrix);
            M0Model temperedModel(&settings, &aln, &temperedTree, &temperedRateMatrix, executor);

            std::vector<Move> temperedMoves;
            temperedMoves.emplace_back(Move{
                settings.treeWeight,
                aln.getNumTaxa(),
                [&temperedTree]() {return temperedTree.update();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.omegaWeight,
                5,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateOmega();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.kWeight,
                3,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateK();},
                []() {return true;}
            });
            temperedMoves.emplace_back(Move{
                settings.stationaryWeight,
                10,
                [&temperedRateMatrix]() {return temperedRateMatrix.updateStationary();},
                []() {return true;}
            });

            MCMC mcmc(&model, moves, &temperedModel, temperedMoves, settings, disableBayesOpt);
            std::cout << "Starting MCMCMC with beta " << settings.mcmcmcBeta
                      << " and swap frequency " << settings.mcmcmcSwapFrequency << "..." << std::endl;
            mcmc.burnin();
            mcmc.run();
        }
        else{
            MCMC mcmc(&model, moves, settings, disableBayesOpt);
            
            std::cout << "Starting MCMC..." << std::endl;
            mcmc.burnin();
            mcmc.run();
        }
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

int main(int argc, char* argv[]){

    Settings settings(argc, argv);
    tf::Executor executor(settings.threads);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    double alphaLength = (settings.treeLengthMean * settings.treeLengthMean) / (settings.treeLengthSD * settings.treeLengthSD);
    double betaLength = settings.treeLengthMean / (settings.treeLengthSD * settings.treeLengthSD);

    double treeLengthParams[2] = {alphaLength, betaLength};

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(settings.simulateM0 == false && settings.simulateCMM == false && settings.simulateDPCMM == false){
        Alignment aln(settings.nexusInput);
        TreeParameter treeParam(aln, settings.tree, treeLengthParams);
        inference(settings, aln, treeParam, settings.bayesOpt == 0, executor);
    }
    else{
        int numSites = 200;
        int numTaxa = 25;
        std::string mcmcOutput = settings.mcmcOutput;
        std::string tipsOutput = settings.tipsOutput;
        std::string dppOutput = settings.dppOutput;
        std::string ancestralStatesOutput = settings.ancestralStatesOutput;
        std::string treeOutput = settings.treeOutput;
        std::string branchOutput = settings.branchOutput;
        int fixedRegimes = settings.fixedRegimes;
        std::vector<int> fixedDPAssignments = settings.fixedDPAssignments;
        
        std::string modelClass = "";
        if(settings.simulateM0)
            modelClass = "M0";
        else if(settings.simulateCMM)
            modelClass = "CMM";
        else if(settings.simulateDPCMM)
            modelClass = "DPCMM";

        std::cout << "Doing a simulation analysis of " << modelClass << " with " << settings.numSimulations << " simulations of " << numTaxa << " taxa trees with " << numSites << " sites." << std::endl;
        for(int i = 0; i < settings.numSimulations; i++){
            std::cout << "Starting simulation " << i << "..." << std::endl;

            TreeObject tree(numTaxa, true);
            double treeLength = Probability::Gamma::rv(&rng, alphaLength, betaLength);
            tree.setTreeLength(treeLength);
            std::vector<Node*> preOrderSeq = tree.getPostOrderSeq();
            std::reverse(preOrderSeq.begin(), preOrderSeq.end());
            Node* root = tree.getRoot();
            int numNodes = tree.getNumNodes();

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
            int trueActiveOmegas = 0;
            std::vector<int> trueDPAssignments;

            if(settings.simulateM0){
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
                            for(Node* d : n->getNeighborRef()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    const Matrix<double>& P = transProb(0, 0, dIndex);
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
                fs << "TreeLength\tK\tOmega";
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
                fs << treeLength << "\t" << rateMatrix.getK() << "\t" << rateMatrix.getOmega();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << stationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs << std::endl;
                fs.close();
            }
            else if(settings.simulateCMM){
                TransitionProbability transProb(numNodes, 305);
                CMMMatrix rateMatrix(settings);
                int activeOmegas = (int)(rng.uniformRv() * 5) + 1;
                trueActiveOmegas = activeOmegas;

                rateMatrix.setActiveOmegas(activeOmegas);
                std::vector<double> stationary = rateMatrix.getStationary();
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                transProb.updateQ(rateMatrix.Q(), 0);
                std::array<double, 5> dNdS = rateMatrix.dNdS();

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
                            for(Node* d : n->getNeighborRef()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    const Matrix<double>& P = transProb(0, 0, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            int state = sites[c*numNodes + nIndex];
                            truedNdS[c*numTaxa + nIndex] = dNdS[(int)(state/61.0)];
                            tipSites[c*numTaxa + nIndex] = state;
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "TreeLength\tOmegaCount\tK\tOmega\tOmegaIncrement1\tOmegaIncrement2\tOmegaIncrement3\tOmegaIncrement4\tR";
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
                fs << treeLength << "\t" << activeOmegas << "\t" << rateMatrix.getK() << "\t" << rateMatrix.getOmega(0) 
                   << "\t" << rateMatrix.getOmega(1) << "\t" << rateMatrix.getOmega(2) << "\t" 
                   << rateMatrix.getOmega(3) << "\t" << rateMatrix.getOmega(4) << "\t" << rateMatrix.getR();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << rawStationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs << std::endl;
                fs.close();
            }
            else if(settings.simulateDPCMM){
                TransitionProbability transProb(numNodes, 183);
                DPCMMMatrix rateMatrix(settings, numSites);
                int activeOmegas = (int)(rng.uniformRv() * 3) + 1;
                trueActiveOmegas = activeOmegas;

                rateMatrix.setActiveOmegas(activeOmegas);
                std::vector<double> stationary = rateMatrix.getStationary();
                std::vector<double> rawStationary = rateMatrix.getRawStationary();
                std::vector<Category> categories = rateMatrix.getCategories();
                std::vector<int> assignments = rateMatrix.getAssignments();
                trueDPAssignments = assignments;
                std::vector<std::array<double, 3>> dNdSVec;

                transProb.allocateQ(categories.size());

                for(int cat = 0; cat < categories.size(); cat++){
                    transProb.updateQ(rateMatrix.Q(cat), cat);
                    std::array<double, 3> dNdS = rateMatrix.dNdS(cat);
                    dNdSVec.push_back(dNdS);
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
                            for(Node* d : n->getNeighborRef()){
                                if(d != n->getAncestor()){
                                    int dIndex = d->getIndex();
                                    const Matrix<double>& P = transProb(0, cat, dIndex);
                                    sites[c*numNodes + dIndex] = randomTransition(rng, ancestralState, P);
                                }
                            }
                        }
                        else{
                            int state = sites[c*numNodes + nIndex];
                            truedNdS[c*numTaxa + nIndex] = dNdSVec[cat][(int)(state/61.0)];
                            tipSites[c*numTaxa + nIndex] = state;
                        }
                    }
                }
                std::cout << "Writing true parameters to file..." << std::endl;
                std::ofstream fs;
                fs.open(settings.simulationOutput + std::to_string(i), std::ofstream::out);
                fs << "TreeLength\tOmegaCount\tK\tR\tCategoryCount";
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
                fs << treeLength << "\t" << activeOmegas << "\t" << rateMatrix.getK() << "\t" << rateMatrix.getR() << "\t" << categories.size();
                for(int element = 0; element < 61; element++){
                    fs << "\t" << rawStationary[element];
                }
                fs << "\t" << tree.getNewick();
                for(int taxon = 0; taxon < numTaxa; taxon++){
                    for(int site = 0; site < numSites; site++){
                        fs << "\t" << truedNdS[site*numTaxa + taxon];
                    }
                }
                fs << std::endl;
                fs.close();
            }
        
            std::cout << "Simulation " << i << " is complete. Starting inference..." << std::endl;
            Alignment aln(tipSites, numSites, numTaxa);
            delete [] sites;
            delete [] tipSites;
            delete [] truedNdS;

            settings.fixedRegimes = fixedRegimes;
            settings.fixedDPAssignments = fixedDPAssignments;
            if(settings.fixCorrectRegime){
                settings.fixedRegimes = trueActiveOmegas;
            }
            if(settings.fixCorrectDP && settings.simulateDPCMM){
                settings.fixedDPAssignments = trueDPAssignments; //Theoretically, we could have users able to set this themselves, but there is no point
            }

            //Change the name I am logging to according to the number simulation it is
            if(!settings.sequentialTuningSim){
                TreeParameter treeParam(tree, treeLengthParams);
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
                
                inference(settings, aln, treeParam, true, executor);
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
                
                TreeParameter treeParamA(tree, treeLengthParams);
                inference(settings, aln, treeParamA, false, executor);

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
                
                TreeParameter treeParamB(tree, treeLengthParams);
                inference(settings, aln, treeParamB, true, executor);
            }
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
}
