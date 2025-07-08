#include "Settings.hpp"
#include "Msg.hpp"
#include <iostream>
#include <string>
#include <vector>

Settings::Settings(int argc,  char* argv[]) : nexusInput(""), treeOutput(""), dppOutput(""), mcmcOutput(""), tipsOutput(""),
                                              branchOutput(""), numIterations(10000), printFrequency(10), sampleFrequency(25),
                                              burnInIterations(1000), tuneFrequency(100), rLambda(5.0), gammaLambda(5.0),
                                              kLambda(2.0), omegaLambda(2.0), treeLengthLambda(5.0), expectedCat(1.2),
                                              rWeight(1.0), kWeight(1.0), stationaryWeight(2.0), omegaWeight(1.0), proportionsWeight(1.0),
                                              dppWeight(2.0), treeWeight(2.0), rjWeight(1.0), numGibbs(5), fixedTree(""), M0(false),
                                              M3S2(false), simulateDPP(false), simulateM0(false), simulateM3S2(false),
                                              numSimulations(1), bayesOpt(0), bayesOptFrequency(0), sequentialTuningSim(false),
                                              SB(false), simulateSB(false), RJ(false), RJDPP(false), truncation(5) {

    std::vector<std::string> settings;
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i];
        settings.push_back(arg);
    }
    #if TEST_RUN==1
    settings.push_back("-treeOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/res/trees.trees");
    settings.push_back("-mcmcOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/res/analysis.log");
    settings.push_back("-dppOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/res/dpp.log");
    settings.push_back("-tipsOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/res/tips.log");
    settings.push_back("-branchOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/res/branches.log");
    settings.push_back("-simulationOutput");
    settings.push_back("/workspaces/Varying_Selection_DPP/res/simulation.log");
    settings.push_back("-simulateM3S2");
    settings.push_back("-RJDPP");
    #endif
    if (settings.size() == 0) {
        usage();
        Msg::error("Expected command line arguments");
    }

    std::string currentArg = "";
    for (int i=0; i<settings.size(); i++) {
        if(settings[i] == "-M0"){
            M0 = true;
            if(M3S2 || SB || RJ || RJDPP){
                Msg::error("Cannot do inference under two models!");
            }
        }
        else if(settings[i] == "-M3S2"){
            M3S2 = true;
            if(M0 || SB || RJ || RJDPP){
                Msg::error("Cannot do inference under two models!");
            }
        }
        else if(settings[i] == "-SB"){
            SB = true;
            Msg::warning("The stick breaking model is experimental! This could go badly!");
            if(M0 || M3S2 || RJ || RJDPP){
                Msg::error("Cannot do inference under two models!");
            }
        }
        else if(settings[i] == "-RJ"){
            RJ = true;
            if(M0 || M3S2 || SB || RJDPP){
                Msg::error("Cannot do inference under two models!");
            }
        }
        else if(settings[i] == "-RJDPP"){
            RJDPP = true;
            if(M0 || M3S2 || SB || RJ){
                Msg::error("Cannot do inference under two models!");
            }
        }
        else if(settings[i] == "-simulateM0"){
            simulateM0 = true;
            if(simulateM3S2 || simulateDPP || simulateSB){
                Msg::error("Cannot simulate under multiple models!");
            }
        }
        else if(settings[i] == "-simulateM3S2"){
            simulateM3S2 = true;
            if(simulateM0 || simulateDPP || simulateSB){
                Msg::error("Cannot simulate under multiple models!");
            }
        }
        else if(settings[i] == "-simulateDPP"){
            simulateDPP = true;
            if(simulateM3S2 || simulateM0 || simulateSB){
                Msg::error("Cannot simulate under multiple models!");
            }
        }
        else if(settings[i] == "-simulateSB"){
            simulateSB = true;
            if(simulateM3S2 || simulateM0 || simulateDPP){
                Msg::error("Cannot simulate under multiple models!");
            }
        }
        else if(settings[i] == "-sequentialTuningSim"){
            sequentialTuningSim = true;
        }
        else if (currentArg == "")
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
            else if (currentArg == "-tipsOut")
                tipsOutput = settings[i];
            else if (currentArg == "-branchOut")
                branchOutput = settings[i];
            else if (currentArg == "-ancestralStatesOut")
                ancestralStatesOutput = settings[i];
            else if (currentArg == "-numIter")
                numIterations = stoi(settings[i]);
            else if (currentArg == "-numGibbs")
                numGibbs = stoi(settings[i]);
            else if (currentArg == "-printFreq")
                printFrequency = stoi(settings[i]);
            else if (currentArg == "-sampleFreq")
                sampleFrequency = stoi(settings[i]);
            else if (currentArg == "-burnInIter")
                burnInIterations = stoi(settings[i]);
            else if (currentArg == "-tuneFreq")
                tuneFrequency = stoi(settings[i]);
            else if (currentArg == "-bayesOpt")
                bayesOpt = stoi(settings[i]);
            else if (currentArg == "-bayesOptFreq")
                bayesOptFrequency = stoi(settings[i]);
            else if (currentArg == "-treeLambda")
                treeLengthLambda = stod(settings[i]);
            else if (currentArg == "-omegaLambda")
                omegaLambda = stod(settings[i]);
            else if (currentArg == "-kLambda")
                kLambda = stod(settings[i]);
            else if (currentArg == "-rLambda")
                rLambda = stod(settings[i]);
            else if (currentArg == "-gammaLambda")
                gammaLambda = stod(settings[i]);
            else if (currentArg == "-expectedCat")
                expectedCat = stod(settings[i]);
             else if (currentArg == "-dppWeight")
                dppWeight = stod(settings[i]);
            else if (currentArg == "-kWeight")
                kWeight = stod(settings[i]);
            else if (currentArg == "-rWeight")
                rWeight = stod(settings[i]);
            else if (currentArg == "-stationaryWeight")
                stationaryWeight = stod(settings[i]);
            else if (currentArg == "-treeWeight")
                treeWeight = stod(settings[i]);
            else if (currentArg == "-omegaWeight")
                omegaWeight = stod(settings[i]);
            else if (currentArg == "-proportionsWeight")
                proportionsWeight = stod(settings[i]);
            else if (currentArg == "-stickBreakingTruncation")
                truncation = stoi(settings[i]);
            else if (currentArg == "-fixedTree")
                fixedTree = settings[i];
            else if (currentArg == "-numSimulations")
                numSimulations = stoi(settings[i]);
            else if (currentArg == "-simulationOutput")
                simulationOutput = settings[i];
            else{
                Msg::error("Could not interpret argument " + currentArg);
                usage();
            }
            currentArg = "";
        }
    }

    bool simulating = (simulateDPP == true || simulateM0 == true || simulateM3S2 == true || simulateSB == true);

    if((nexusInput == "" || treeOutput == "" || mcmcOutput == "") && !simulating){
        usage();
        Msg::error("For non-simulation analyses, nexusInput, treeOut, mcmcOut are required arguments.");
    }

    if(simulating && nexusInput != ""){
        usage();
        Msg::warning("Simulation analyses cannot use a nexus input. This file will be ignored!");
    }

    if(simulating && fixedTree != ""){
        usage();
        Msg::warning("Simulation analyses cannot use a provided tree. This file will be ignored!");
    }

    if(simulating && simulationOutput == ""){
        usage();
        Msg::warning("Simulation analyses require an output to be set!");
    }

    if(simulating && (treeOutput == "" || mcmcOutput == "")){
        Msg::error("For simulation analyses, treeOut, mcmcOut are required arguments.");
    }

    if(sequentialTuningSim && !simulating){
        Msg::error("For a sequential tuning simulation, the software must be in simulation mode.");
    }

    print();
}


void Settings::print(){
    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : " << nexusInput << std::endl;
    std::cout << "   * -treeOut           : " << treeOutput << std::endl;
    std::cout << "   * -mcmcOut           : " << mcmcOutput << std::endl;
    std::cout << "   * -dppOut            : " << dppOutput << std::endl;
    std::cout << "   * -tipsOut           : " << tipsOutput << std::endl;
    std::cout << "   * -ancestralStatesOut: " << ancestralStatesOutput << std::endl;
    std::cout << "   * -fixedTree         : " << fixedTree << std::endl;
    std::cout << "   * -simulationOutput  : " << simulationOutput << std::endl;
    std::cout << std::endl;

    std::cout << "Inference Model and Simulation:" << std::endl;
    std::cout << "   * -M0                : " << M0 << std::endl;
    std::cout << "   * -M3S2              : " << M3S2 << std::endl;
    std::cout << "   * -SB                : " << SB << std::endl;
    std::cout << "   * -RJ                : " << RJ << std::endl;
    std::cout << "   * -RJDPP             : " << RJDPP << std::endl;
    std::cout << "   * -simulateM0        : " << simulateM0 << std::endl;
    std::cout << "   * -simulateM3S2      : " << simulateM3S2 << std::endl;
    std::cout << "   * -simulateDPP       : " << simulateDPP << std::endl;
    std::cout << "   * -simulateSB        : " << simulateSB << std::endl;
    std::cout << "   * -numSimulations    : " << numSimulations << std::endl;
    std::cout << std::endl;
    
    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : " << treeLengthLambda << std::endl;
    std::cout << "   * -omegaLambda       : " << omegaLambda << std::endl;
    std::cout << "   * -kLambda           : " << kLambda << std::endl;
    std::cout << "   * -rLambda           : " << rLambda << std::endl;
    std::cout << "   * -gammaLambda       : " << gammaLambda << std::endl;
    std::cout << "   * -expectedCat       : " << expectedCat << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : " << numIterations << std::endl;
    std::cout << "   * -printFreq         : " << printFrequency << std::endl;
    std::cout << "   * -sampleFreq        : " << sampleFrequency << std::endl;
    std::cout << "   * -burnInIter        : " << burnInIterations << std::endl;
    std::cout << "   * -tuneFreq          : " << tuneFrequency << std::endl;
    std::cout << "   * -bayesOpt          : " << bayesOpt << std::endl;
    std::cout << "   * -bayesOptFreq      : " << bayesOptFrequency << std::endl;
    std::cout << "   * -numGibbs          : " << numGibbs << std::endl;
    std::cout << "   * -treeWeight        : " << treeWeight << std::endl;
    std::cout << "   * -kWeight           : " << kWeight << std::endl;
    std::cout << "   * -rWeight           : " << rWeight << std::endl;
    std::cout << "   * -rjWeight          : " << rjWeight << std::endl;
    std::cout << "   * -stationaryWeight  : " << stationaryWeight << std::endl;
    std::cout << "   * -dppWeight         : " << dppWeight << std::endl;
    std::cout << "   * -omegaWeight       : " << omegaWeight << std::endl;
    std::cout << std::endl;
}

void Settings::usage(void) {

    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : Input nexus file containing the nculeotide alignment." << std::endl;
    std::cout << "   * -treeOut           : The output file name for the tree trace." << std::endl;
    std::cout << "   * -mcmcOut           : The output file name for the bulk of the MCMC trace, excluding the tree and DPP parameters." << std::endl;
    std::cout << "   * -dppOut            : The output file name for the DPP parameters." << std::endl;
    std::cout << "   * -tipsOut           : The output file name for the reconstructed tip dNdS ratios." << std::endl;
    std::cout << "   * -ancestralStatesOut: The output file name for the all ancestral dNdS ratios." << std::endl;
    std::cout << "   * -fixedTree         : The NEWICK string corresponding to the fixed tree you wish to analyze." << std::endl;
    std::cout << "   * -simulationOutput  : The output file name for the true simulation parameters." << std::endl;
    std::cout << std::endl;

    std::cout << "Inference Model and Simulation:" << std::endl;
    std::cout << "   * -M0                : Do inference under a normal codon phylogenetic model." << std::endl;
    std::cout << "   * -M3S2              : Do inference under M3S2 as described by Guindon et al. (2004)." << std::endl;
    std::cout << "   * -SB                : (Experimental) Do inference under a stick-breaking Markov-modulated model." << std::endl;
    std::cout << "   * -RJ                : (Experimental) Do inference under the reversible jump Markov-modulated model." << std::endl;
    std::cout << "   * -RJDPP             : (Experimental) Do inference with the reversible-jump DPP model." << std::endl;
    std::cout << "   * -simulateM0        : Directs the program to simulate under M0 and test against the selected inference model." << std::endl;
    std::cout << "   * -simulateM3S2      : Directs the program to simulate under M3S2 and test against the selected inference model." << std::endl;
    std::cout << "   * -simulateDPP       : Directs the program to simulate under the DPP model and test against the selected inference model." << std::endl;
    std::cout << "   * -simulateSB        : (Experimental) Directs the program to simulate under the stick-breaking Markov-modulated model and test against the selected inference model." << std::endl;
    std::cout << "   * -numSimulations    : The number of simulations to do inference under." << std::endl;
    std::cout << std::endl;

    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : Lambda parameter for the tree length exponential prior." << std::endl;
    std::cout << "   * -omegaLambda       : Lambda parameter for the nonsynonymous mutation rate's exponential prior." << std::endl;
    std::cout << "   * -kLambda           : Lambda parameter for the transition/transversion rate's exponential prior." << std::endl;
    std::cout << "   * -rLambda           : Lambda parameter for the matrix-swapping rate's exponential prior." << std::endl;
    std::cout << "   * -gammaLambda       : Lambda parameter for the global matrix-swapping rate's exponential prior." << std::endl;
    std::cout << "   * -expectedCat       : The number of expected categories for the DPP." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : The number of iterations for the MCMC." << std::endl;
    std::cout << "   * -printFreq         : How often to output the MCMC state to the screen." << std::endl;
    std::cout << "   * -sampleFreq        : How often to ouput the MCMC state to log files." << std::endl;
    std::cout << "   * -burnInIter        : The number of iterations for the burn-in." << std::endl;
    std::cout << "   * -tuneFreq          : How often to tune the acceptance rate of the MCMC moves during the burn-in." << std::endl;
    std::cout << "   * -bayesOpt          : (Experimental) The number of iterations to run Bayesian optimization on the MCMC moves after the burn-in." << std::endl;
    std::cout << "   * -bayesOptFreq      : (Experimental) How often to sample the trace for Bayesian optimization." << std::endl;
    std::cout << "   * -numGibbs          : How many Gibbs updates to perform on the DPP partitions." << std::endl;
    std::cout << "   * -treeWeight        : How often to propose a move on the tree." << std::endl;
    std::cout << "   * -kWeight           : How often to propose a move on the K parameter." << std::endl;
    std::cout << "   * -rWeight           : How often to propose a move on the R or Gamma parameter." << std::endl;
    std::cout << "   * -rWeight           : How often to propose a move on the dimension of the Markov modulated model." << std::endl;
    std::cout << "   * -stationaryWeight  : How often to propose a move on the stationary distribution." << std::endl;
    std::cout << "   * -dppWeight         : How often to propose a move on the DPP partitions." << std::endl;
    std::cout << "   * -omegaWeight       : How often to propose a move on the omega parameters." << std::endl;
    std::cout << std::endl;
}