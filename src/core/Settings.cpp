#include "Settings.hpp"
#include "Msg.hpp"
#include <iostream>
#include <string>
#include <vector>

Settings::Settings(int argc,  char* argv[]) : nexusInput(""), treeOutput(""), dppOutput(""), mcmcOutput(""),
                                              numIterations(100000), printFrequency(10), sampleFrequency(500),
                                              burnInIterations(10000), tuneFrequency(2500), rLambda(2.0),
                                              kLambda(1.0), omegaLambda(1.0), dppAlpha(0.5), updateStationary(true),
                                              numGibbsUpdate(10), rateMatrixWeight(10), dppWeight(10), treeWeight(5),
                                              treeLengthLambda(1.0), simulate(false), fixedTree(""), numTaxa(-1),
                                              numChar(-1), kValue(-1.0), rValue(-1.0), scaleTree(1.0) {
    std::vector<std::string> settings;
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i];
        settings.push_back(arg);
    }

    /*
    settings.push_back("-nexus");
    settings.push_back("/workspaces/Varying_Selection_DPP/publication_analyses/globin_analysis/globins.nex");
    settings.push_back("-treeOut");
    settings.push_back("/workspaces/Code/Varying_Selection_DPP/res/trees.trees");
    settings.push_back("-mcmcOut"); 
    settings.push_back("/workspaces/Code/Varying_Selection_DPP/res/analysis.log");
    settings.push_back("-dppOut"); 
    settings.push_back("/workspaces/Code/Varying_Selection_DPP/res/dpp.log");
    settings.push_back("-fixedTree");
    settings.push_back("((((CcarpioHBB:0.05172908,DrerioHBB:0.03752478):0.07009799,SsalarHBB:0.14769959):0.19916559,XborealisHBB:0.47154023):0.06222666,(GgallusHBBA:0.13597721,(BtaurusHBB:0.08606479,HsapiensHBB:0.06002410):0.13346986):0.05790206,(((CcarpioHBA:0.09708301,DrerioHBA:0.02472006):0.15923891,SsalarHBA:0.20586909):0.16195969,(XborealisHBA:0.38418109,(GgallusHBAA:0.12817884,(BtaurusHBA2:0.06468952,HsapiensHBA2:0.07819764):0.12209754):0.07038107):0.06943181):0.22254004);");
    */

/*
    settings.push_back("-simulate");
    settings.push_back("1");
    settings.push_back("-omega1Vector");
    settings.push_back("2.0,0.5");
    settings.push_back("-numChar");
    settings.push_back("250");
    settings.push_back("-omega2Vector");
    settings.push_back("1.0,0.25");
    settings.push_back("-kValue");
    settings.push_back("2.5");
    settings.push_back("-rValue");
    settings.push_back("0.15");
    settings.push_back("-treeLambda");
    settings.push_back("0.5");
*/

    if (settings.size() == 0) {
        usage();
        Msg::error("Expected command line arguments");
    }

    std::string currentArg = "";
    for (int i=0; i<settings.size(); i++) {
        if (currentArg == "")
            currentArg = settings[i];
        else {
            if (currentArg == "-nexus")
                nexusInput = settings[i];
            else if (currentArg == "-treeOut")
                treeOutput = settings[i];
            else if (currentArg == "-mcmcOut")
                mcmcOutput = settings[i];
            else if (currentArg == "-dppOut")
                dppOutput = settings[i];
            else if (currentArg == "-numIter")
                numIterations = stoi(settings[i]);
            else if (currentArg == "-printFreq")
                printFrequency = stoi(settings[i]);
            else if (currentArg == "-sampleFreq")
                sampleFrequency = stoi(settings[i]);
            else if (currentArg == "-burnInIter")
                burnInIterations = stoi(settings[i]);
            else if (currentArg == "-tuneFreq")
                tuneFrequency = stoi(settings[i]);
            else if (currentArg == "-treeLambda")
                treeLengthLambda = stod(settings[i]);
            else if (currentArg == "-omegaLambda")
                omegaLambda = stod(settings[i]);
            else if (currentArg == "-kLambda")
                kLambda = stod(settings[i]);
            else if (currentArg == "-rLambda")
                rLambda = stod(settings[i]);
            else if (currentArg == "-dppAlpha")
                dppAlpha = stod(settings[i]);
             else if (currentArg == "-dppWeight")
                dppWeight = stod(settings[i]);
            else if (currentArg == "-rateMatrixWeight")
                rateMatrixWeight = stod(settings[i]);
            else if (currentArg == "-treeWeight")
                treeWeight = stod(settings[i]);
            else if (currentArg == "-updateStationary")
                updateStationary = stoi(settings[i]) == 1;
            else if (currentArg == "-numGibbsUpdate")
                numGibbsUpdate = stoi(settings[i]);
            else if (currentArg == "-simulate")
                simulate = stoi(settings[i]) == 1;
            else if (currentArg == "-fixedTree")
                fixedTree = settings[i];
            else if (currentArg == "-numTaxa")
                numTaxa = stoi(settings[i]);
            else if (currentArg == "-numChar")
                numChar = stoi(settings[i]);
            else if (currentArg == "-kValue")
                kValue = stod(settings[i]);
            else if (currentArg == "-rValue")
                rValue = stod(settings[i]);
            else if (currentArg == "-omega1Vector"){
                std::string currentString = "";
                for(int c = 0; c < settings[i].size(); c++){
                    if(settings[i][c] == ','){
                        omega1Vector.push_back(stod(currentString));
                        currentString = "";
                    }
                    else{
                        currentString += settings[i][c];
                    }
                }
                omega1Vector.push_back(stod(currentString));
            }
            else if (currentArg == "-omega2Vector"){
                std::string currentString = "";
                for(int c = 0; c < settings[i].size(); c++){
                    if(settings[i][c] == ','){
                        omega2Vector.push_back(stod(currentString));
                        currentString = "";
                    }
                    else{
                        currentString += settings[i][c];
                    }
                }
                omega2Vector.push_back(stod(currentString));
            }
            else if (currentArg == "-assignmentVector"){
                std::string currentString = "";
                for(int c = 0; c < settings[i].size(); c++){
                    if(settings[i][c] == ','){
                        assignmentVector.push_back(stoi(currentString));
                        currentString = "";
                    }
                    else{
                        currentString += settings[i][c];
                    }
                }
                assignmentVector.push_back(stoi(currentString));
            }
            else{
                Msg::error("Could not interpret argument " + currentArg);
                usage();
            }
            currentArg = "";
        }
    }

    if(!simulate){
        if(nexusInput == "" || treeOutput == "" || mcmcOutput == "" || dppOutput == ""){
            usage();
            Msg::error("For non-simulation analyses, nexusInput, treeOutput, mcmcOutput, and dppOutput are required arguments.");
        }
    }
    else{
        if(omega1Vector.size() == 0 || omega2Vector.size() == 0){
            usage();
            Msg::error("Simulation analyses require the omega vectors to be defined.");
        }
    }

    print();
}


void Settings::print(){
    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : " << nexusInput << std::endl;
    std::cout << "   * -treeOut           : " << treeOutput << std::endl;
    std::cout << "   * -mcmcOut           : " << mcmcOutput << std::endl;
    std::cout << "   * -dppOut            : " << dppOutput << std::endl;
    std::cout << "   * -fixedTree         : " << fixedTree << std::endl;
    std::cout << std::endl;
    
    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : " << treeLengthLambda << std::endl;
    std::cout << "   * -omegaLambda       : " << omegaLambda << std::endl;
    std::cout << "   * -kLamdba           : " << kLambda << std::endl;
    std::cout << "   * -rLambda           : " << rLambda << std::endl;
    std::cout << "   * -dppAlpha          : " << dppAlpha << std::endl;
    std::cout << "   * -kValue            : " << kValue << std::endl;
    std::cout << "   * -rValue            : " << rValue << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : " << numIterations << std::endl;
    std::cout << "   * -printFreq         : " << printFrequency << std::endl;
    std::cout << "   * -sampleFreq        : " << sampleFrequency << std::endl;
    std::cout << "   * -burnInIter        : " << burnInIterations << std::endl;
    std::cout << "   * -tuneFreq          : " << tuneFrequency << std::endl;
    std::cout << "   * -numGibbsUpdate    : " << numGibbsUpdate << std::endl;
    std::cout << "   * -updateStationary  : " << updateStationary << std::endl;
    std::cout << "   * -treeWeight        : " << treeWeight << std::endl;
    std::cout << "   * -rateMatrixWeight  : " << rateMatrixWeight << std::endl;
    std::cout << "   * -dppWeight         : " << dppWeight << std::endl;
    std::cout << std::endl;

    std::cout << "Simulation:" << std::endl;
    std::cout << "   * -simulate          : " << simulate << std::endl;
    std::cout << "   * -numTaxa           : " << numTaxa << std::endl;
    std::cout << "   * -numChar           : " << numChar << std::endl;
    std::cout << "   * -omega1Vector      : ";
    for(double o : omega1Vector)
        std::cout << o << " ";
    std::cout << std::endl;
    std::cout << "   * -omega2Vector      : ";
    for(double o : omega2Vector)
        std::cout << o << " ";
    std::cout << std::endl;
    std::cout << "   * -asignmentVector   : ";
    for(int a : assignmentVector)
        std::cout << a << " ";
    std::cout << std::endl;
}

void Settings::usage(void) {

    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : Input nexus file containing the nculeotide alignment." << std::endl;
    std::cout << "   * -treeOut           : The output file name for the tree trace." << std::endl;
    std::cout << "   * -mcmcOut           : The output file name for the bulk of the MCMC trace, excluding the tree and DPP parameters." << std::endl;
    std::cout << "   * -dppOut            : The output file name for the DPP parameters." << std::endl;
    std::cout << "   * -fixedTree         : The NEWICK string corresponding to the fixed tree you wish to analyze." << std::endl;
    std::cout << std::endl;

    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : Lambda parameter for the tree length exponential prior." << std::endl;
    std::cout << "   * -omegaLambda       : Lambda parameter for the dN/dS exponential prior." << std::endl;
    std::cout << "   * -kLamdba           : Lambda parameter for the transition/transversion ratio exponential prior." << std::endl;
    std::cout << "   * -rLambda           : Lambda parameter for the matrix-swapping exponential prior." << std::endl;
    std::cout << "   * -dppAlpha          : The alpha parameter of the DPP." << std::endl;
    std::cout << "   * -kValue            : The starting value for the K parameter." << std::endl;
    std::cout << "   * -rValue            : The starting value for the R parameter." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : The number of iterations for the MCMC." << std::endl;
    std::cout << "   * -printFreq         : How often to output the MCMC state to the screen." << std::endl;
    std::cout << "   * -sampleFreq        : How often to ouput the MCMC state to log files." << std::endl;
    std::cout << "   * -burnInIter        : The number of iterations for the burn-in." << std::endl;
    std::cout << "   * -tuneFreq          : How often to tune the MCMC moves during the burn-in." << std::endl;
    std::cout << "   * -numGibbsUpdate    : How many sites to use during the Gibbs sampling of the DPP." << std::endl;
    std::cout << "   * -updateStationary  : Whether or not the stationary should be treated as empirical or updated as a parameter (0/1)." << std::endl;
    std::cout << "   * -treeWeight        : How often to propose a move on the tree." << std::endl;
    std::cout << "   * -rateMatrixWeight  : How often to propose a move on the stationary, R, and K." << std::endl;
    std::cout << "   * -dppWeight         : How often to propose a move on the DPP parameters." << std::endl;
    std::cout << std::endl;

    std::cout << "Simulation:" << std::endl;
    std::cout << "   * -simulate          : Should this run be a simulation under the model (0/1)?" << std::endl;
    std::cout << "   * -numTaxa           : The number of taxa to simulate." << std::endl;
    std::cout << "   * -numChar           : The number of chararacters to simulate." << std::endl;
    std::cout << "   * -omega1Vector      : The vector of omega 1 rate multipliers." << std::endl;
    std::cout << "   * -omega2Vector      : The vector of omega 2 rate multipliers." << std::endl;
    std::cout << "   * -asignmentVector   : The vector of category assignments." << std::endl;
    std::cout << std::endl;
}