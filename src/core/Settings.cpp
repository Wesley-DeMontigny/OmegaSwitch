#include "Settings.hpp"
#include "Msg.hpp"
#include <iostream>
#include <string>
#include <vector>

Settings::Settings(int argc,  char* argv[]) : nexusInput(""), treeOutput(""), dppOutput(""), mcmcOutput(""),
                                              numIterations(100000), printFrequency(10), sampleFrequency(100),
                                              burnInIterations(10000), tuneFrequency(500), rLambda(2.0),
                                              kLambda(1.0), omegaLambda(1.0), dppAlpha(0.5), updateStationary(false),
                                              numGibbsUpdate(10), rateMatrixWeight(10), dppWeight(5), treeWeight(10),
                                              treeLengthLambda(10.0) {
    std::vector<std::string> settings;
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i];
        settings.push_back(arg);
    }

    settings.push_back("-nexus");
    settings.push_back("C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/replicase.nex");
    settings.push_back("-treeOut");
    settings.push_back("C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/trees.trees");
    settings.push_back("-mcmcOut"); 
    settings.push_back("C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/analysis.log");
    settings.push_back("-dppOut"); 
    settings.push_back("C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/dpp.log");

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
            else if (currentArg == "-treeLamdba")
                treeLengthLambda = stod(settings[i]);
            else if (currentArg == "-omegaLamdba")
                omegaLambda = stod(settings[i]);
            else if (currentArg == "-kLamdba")
                kLambda = stod(settings[i]);
            else if (currentArg == "-rLamdba")
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
            else
                Msg::error("Could not interpret argument " + settings[i]);
            currentArg = "";
        }
    }

    if(nexusInput == "" || treeOutput == "" || mcmcOutput == "" || dppOutput == ""){
        usage();
        Msg::error("Missing a required argument");
    }

    print();
}


void Settings::print(){
    std::cout << "Input/Output:" << std::endl;
    std::cout << "   * -nexus             : " << nexusInput << std::endl;
    std::cout << "   * -treeOut           : " << treeOutput << std::endl;
    std::cout << "   * -mcmcOut           : " << mcmcOutput << std::endl;
    std::cout << "   * -dppOut            : " << dppOutput << std::endl;
    std::cout << std::endl;
    
    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : " << treeLengthLambda << std::endl;
    std::cout << "   * -omegaLambda       : " << omegaLambda << std::endl;
    std::cout << "   * -kLamdba           : " << kLambda << std::endl;
    std::cout << "   * -rLambda           : " << rLambda << std::endl;
    std::cout << "   * -dppAlpha          : " << dppAlpha << std::endl;
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
}

void Settings::usage(void) {

    std::cout << "Input/Output (Required):" << std::endl;
    std::cout << "   * -nexus             : Input nexus file containing the nculeotide alignment." << std::endl;
    std::cout << "   * -treeOut           : The output file name for the tree trace." << std::endl;
    std::cout << "   * -mcmcOut           : The output file name for the bulk of the MCMC trace, excluding the tree and DPP parameters." << std::endl;
    std::cout << "   * -dppOut            : The output file name for the DPP parameters." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Model Parameters (Optional):" << std::endl;
    std::cout << "   * -treeLambda        : Lambda parameter for the tree length exponential prior." << std::endl;
    std::cout << "   * -omegaLambda       : Lambda parameter for the dN/dS exponential prior." << std::endl;
    std::cout << "   * -kLamdba           : Lambda parameter for the transition/transversion ratio exponential prior." << std::endl;
    std::cout << "   * -rLambda           : Lambda parameter for the matrix-swapping exponential prior." << std::endl;
    std::cout << "   * -dppAlpha          : The alpha parameter of the DPP." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options (Optional):" << std::endl;
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
}