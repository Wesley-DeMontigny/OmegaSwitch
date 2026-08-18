#include "Settings.hpp"
#include "Msg.hpp"
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Settings constructor
 * 
 * @param argc The number of arguments passed from the main function
 * @param argv The arguments passed from the main function
 */
Settings::Settings(int argc,  char* argv[]) {

    std::vector<std::string> settings;
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i];
        settings.push_back(arg);
    }
    #if TEST_RUN==1
    settings.push_back("-mcmcOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/testing/analysis.log");
    settings.push_back("-tipsOut");
    settings.push_back("/workspaces/Varying_Selection_DPP/testing/tips.log");
    settings.push_back("-CMM");
    treeFile = "./testing/primate.treefile";
    nexusInput = "./testing/primate_lys.nex";
    numIterations = 5000;
    burnInIterations = 500;
    tuneFrequency = 100;
    threads = 4;
    #endif
    if(settings.size() == 0) {
        usage();
        Msg::error("Expected command line arguments");
    }

    bool simulating = false;
    bool modelSelected = false;

    std::string currentArg = "";
    for (int i=0; i<settings.size(); i++) {
        if(settings[i] == "-M0"){
            M0 = true;
            if(modelSelected){
                Msg::error("Cannot do inference under two models!");
            }
            modelSelected = true;
        }
        else if(settings[i] == "-CMM"){
            CMM = true;
            if(modelSelected){
                Msg::error("Cannot do inference under two models!");
            }
            modelSelected = true;
        }
        else if(settings[i] == "-DPCMM"){
            DPCMM = true;
            if(modelSelected){
                Msg::error("Cannot do inference under two models!");
            }
            modelSelected = true;
        }
        else if(settings[i] == "-simulateM0"){
            simulateM0 = true;
            if(simulating){
                Msg::error("Cannot simulate under multiple models!");
            }
            simulating = true;
        }
        else if(settings[i] == "-simulateCMM"){
            simulateCMM = true;
            if(simulating){
                Msg::error("Cannot simulate under multiple models!");
            }
            simulating = true;
        }
        else if(settings[i] == "-simulateDPCMM"){
            simulateDPCMM = true;
            if(simulating){
                Msg::error("Cannot simulate under multiple models!");
            }
            simulating = true;
        }
        else if(settings[i] == "-sequentialTuningSim"){
            sequentialTuningSim = true;
        }
        else if(settings[i] == "-fixCorrectRegime"){
            fixCorrectRegime = true;
        }
        else if(settings[i] == "-fixCorrectDP"){
            fixCorrectDP = true;
        }
        else if(currentArg == "")
            currentArg = settings[i];
        else {
            if(currentArg == "-nexus")
                nexusInput = settings[i];
            else if(currentArg == "-treeOut")
                treeOutput = settings[i];
            else if(currentArg == "-mcmcOut")
                mcmcOutput = settings[i];
            else if(currentArg == "-dppOut")
                dppOutput = settings[i];
            else if(currentArg == "-tipsOut")
                tipsOutput = settings[i];
            else if(currentArg == "-branchOut")
                branchOutput = settings[i];
            else if(currentArg == "-ancestralStatesOut")
                ancestralStatesOutput = settings[i];
            else if(currentArg == "-numIter")
                numIterations = stoi(settings[i]);
            else if(currentArg == "-numGibbs")
                numGibbs = stoi(settings[i]);
            else if(currentArg == "-printFreq")
                printFrequency = stoi(settings[i]);
            else if(currentArg == "-sampleFreq")
                sampleFrequency = stoi(settings[i]);
            else if(currentArg == "-burnInIter")
                burnInIterations = stoi(settings[i]);
            else if(currentArg == "-fixedRegimes")
                fixedRegimes = stoi(settings[i]);
            else if(currentArg == "-tuneFreq")
                tuneFrequency = stoi(settings[i]);
            else if(currentArg == "-bayesOpt")
                bayesOpt = stoi(settings[i]);
            else if(currentArg == "-bayesOptFreq")
                bayesOptFrequency = stoi(settings[i]);
            else if(currentArg == "-treeMean")
                treeLengthMean = stod(settings[i]);
            else if(currentArg == "-treeSD")
                treeLengthSD = stod(settings[i]);
            else if(currentArg == "-omegaLambda")
                omegaLambda = stod(settings[i]);
            else if(currentArg == "-kLambda")
                kLambda = stod(settings[i]);
            else if(currentArg == "-rLambda")
                rLambda = stod(settings[i]);
            else if(currentArg == "-expectedCat")
                expectedCat = stod(settings[i]);
            else if(currentArg == "-mcmcmcBeta")
                mcmcmcBeta = stod(settings[i]);
            else if(currentArg == "-mcmcmcSwapFreq")
                mcmcmcSwapFrequency = stoi(settings[i]);
             else if(currentArg == "-dppWeight")
                dppWeight = stod(settings[i]);
            else if(currentArg == "-kWeight")
                kWeight = stod(settings[i]);
            else if(currentArg == "-rWeight")
                rWeight = stod(settings[i]);
            else if(currentArg == "-stationaryWeight")
                stationaryWeight = stod(settings[i]);
            else if(currentArg == "-treeWeight")
                treeWeight = stod(settings[i]);
            else if(currentArg == "-omegaWeight")
                omegaWeight = stod(settings[i]);
            else if(currentArg == "-tree")
                tree = settings[i];
            else if(currentArg == "-treeFile")
                treeFile = settings[i];
            else if(currentArg == "-numSimulations")
                numSimulations = stoi(settings[i]);
            else if(currentArg == "-simulationOutput")
                simulationOutput = settings[i];
            else if(currentArg == "-threads")
                threads = stoi(settings[i]);
            else{
                Msg::error("Could not interpret argument " + currentArg);
                usage();
            }
            currentArg = "";
        }
    }

    bool noTree = tree == "" && treeFile == "";

    if((nexusInput == "" || mcmcOutput == "" || noTree) && !simulating){
        usage();
        Msg::error("For non-simulation analyses, nexus, mcmcOut, and a topology are required arguments.");
    }

    if(simulating && nexusInput != ""){
        usage();
        Msg::warning("Simulation analyses cannot use a nexus input. This file will be ignored!");
    }

    if(simulating && !noTree){
        usage();
        Msg::warning("Simulation analyses cannot use a provided tree. This file be ignored!");
    }

    if(simulating && simulationOutput == ""){
        usage();
        Msg::error("Simulation analyses require an output to be set!");
    }

    if(sequentialTuningSim && !simulating){
        Msg::error("For a sequential tuning simulation, the software must be in simulation mode.");
    }

    if(fixCorrectRegime && !simulating){
        Msg::error("-fixCorrectRegime can only be used during simulation analyses.");
    }

    if(fixCorrectRegime && !simulateCMM && !simulateDPCMM){
        Msg::error("-fixCorrectRegime requires simulating under CMM or DPCMM.");
    }

    if(fixCorrectRegime && !CMM && !DPCMM){
        Msg::error("-fixCorrectRegime requires inference under CMM or DPCMM.");
    }

    if(fixCorrectDP && !simulating){
        Msg::error("-fixCorrectDP can only be used during simulation analyses.");
    }

    if(fixCorrectDP && !simulateDPCMM){
        Msg::error("-fixCorrectDP requires simulating under DPCMM.");
    }

    if(fixCorrectDP && !DPCMM){
        Msg::error("-fixCorrectDP requires DPCMM inference.");
    }

    if(simulating && numSimulations == 0)
        numSimulations = 1;

    if(mcmcmcBeta <= 0.0 || mcmcmcBeta > 1.0){
        Msg::error("The MCMCMC beta must be in the interval (0, 1].");
    }

    if(mcmcmcSwapFrequency < 1){
        Msg::error("The MCMCMC swap frequency must be at least 1.");
    }

    if(fixedRegimes < 0){
        Msg::error("The fixed regime count must be non-negative.");
    }

    if(fixedRegimes > 0){
        if(CMM){
            if(fixedRegimes < 1 || fixedRegimes > 5){
                Msg::error("The fixed regime count for CMM must be in the interval [1, 5].");
            }
        }
        else if(DPCMM){
            if(fixedRegimes < 1 || fixedRegimes > 3){
                Msg::error("The fixed regime count for DPCMM must be in the interval [1, 3].");
            }
        }
        else{
            Msg::error("A fixed regime count can only be used with the CMM or DPCMM inference models.");
        }
    }

    print();
}

/**
 * @brief Prints out the parameters values and settings that are being used this run 
 * 
 */
void Settings::print(){
    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : " << nexusInput << std::endl;
    std::cout << "   * -treeOut           : " << treeOutput << std::endl;
    std::cout << "   * -mcmcOut           : " << mcmcOutput << std::endl;
    std::cout << "   * -dppOut            : " << dppOutput << std::endl;
    std::cout << "   * -tipsOut           : " << tipsOutput << std::endl;
    std::cout << "   * -ancestralStatesOut: " << ancestralStatesOutput << std::endl;
    std::cout << "   * -tree              : " << tree << std::endl;
    std::cout << "   * -treeFile          : " << treeFile << std::endl;
    std::cout << "   * -simulationOutput  : " << simulationOutput << std::endl;
    std::cout << "   * -threads           : " << threads << std::endl;
    std::cout << std::endl;

    std::cout << "Inference Model and Simulation:" << std::endl;
    std::cout << "   * -M0                : " << M0 << std::endl;
    std::cout << "   * -CMM               : " << CMM << std::endl;
    std::cout << "   * -DPCMM             : " << DPCMM << std::endl;
    std::cout << "   * -simulateM0        : " << simulateM0 << std::endl;
    std::cout << "   * -simulateCMM       : " << simulateCMM << std::endl;
    std::cout << "   * -simulateDPCMM     : " << simulateDPCMM << std::endl;
    std::cout << "   * -fixCorrectRegime  : " << fixCorrectRegime << std::endl;
    std::cout << "   * -fixCorrectDP      : " << fixCorrectDP << std::endl;
    std::cout << "   * -numSimulations    : " << numSimulations << std::endl;
    std::cout << std::endl;
    
    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeMean          : " << treeLengthMean << std::endl;
    std::cout << "   * -treeSD            : " << treeLengthSD << std::endl;
    std::cout << "   * -omegaLambda       : " << omegaLambda << std::endl;
    std::cout << "   * -kLambda           : " << kLambda << std::endl;
    std::cout << "   * -rLambda           : " << rLambda << std::endl;
    std::cout << "   * -expectedCat       : " << expectedCat << std::endl;
    std::cout << "   * -fixedRegimes      : " << fixedRegimes << std::endl;
    std::cout << "   * -mcmcmcBeta        : " << mcmcmcBeta << std::endl;
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
    std::cout << "   * -mcmcmcSwapFreq    : " << mcmcmcSwapFrequency << std::endl;
    std::cout << "   * -treeWeight        : " << treeWeight << std::endl;
    std::cout << "   * -kWeight           : " << kWeight << std::endl;
    std::cout << "   * -rWeight           : " << rWeight << std::endl;
    std::cout << "   * -rjWeight          : " << rjWeight << std::endl;
    std::cout << "   * -stationaryWeight  : " << stationaryWeight << std::endl;
    std::cout << "   * -dppWeight         : " << dppWeight << std::endl;
    std::cout << "   * -omegaWeight       : " << omegaWeight << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Informs the user how to use our software from the command line
 * 
 */
void Settings::usage() {

    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : Input nexus file containing the nculeotide alignment." << std::endl;
    std::cout << "   * -treeOut           : The output file name for the tree trace." << std::endl;
    std::cout << "   * -mcmcOut           : The output file name for the bulk of the MCMC trace, excluding the tree and DP parameters." << std::endl;
    std::cout << "   * -dppOut            : The output file name for the DP parameters." << std::endl;
    std::cout << "   * -tipsOut           : The output file name for the reconstructed tip dN/dS ratios." << std::endl;
    std::cout << "   * -ancestralStatesOut: The output file name for the all ancestral dN/dS ratios." << std::endl;
    std::cout << "   * -tree              : The NEWICK string corresponding to the fixed tree you wish to analyze." << std::endl;
    std::cout << "   * -treeFile          : The file containing the NEWICK tree string for the fixed tree you wish to analyze." << std::endl;
    std::cout << "   * -simulationOutput  : The output file name for the true simulation parameters." << std::endl;
    std::cout << "   * -threads           : The number of threads to use during the analysis." << std::endl;
    std::cout << std::endl;

    std::cout << "Inference Model and Simulation:" << std::endl;
    std::cout << "   * NOTE: By default, this software will run the M0 model." << std::endl;
    std::cout << "   * -M0                : Do inference under a normal codon phylogenetic model." << std::endl;
    std::cout << "   * -CMM               : Do inference under the reversible jump Markov-modulated model." << std::endl;
    std::cout << "   * -DPCMM             : Do inference with the reversible-jump DP model." << std::endl;
    std::cout << "   * -simulateM0        : Directs the program to simulate under M0 and test against the selected inference model." << std::endl;
    std::cout << "   * -simulateCMM       : Directs the program to simulate under the CMM model and test against the selected inference model." << std::endl;
    std::cout << "   * -simulateDPCMM     : Directs the program to simulate under the DPCMM model and test against the selected inference model." << std::endl;
    std::cout << "   * -fixCorrectRegime  : During CMM or DPCMM simulations, fix inference to the true simulated regime count." << std::endl;
    std::cout << "   * -fixCorrectDP      : During DPCMM simulations, fix inference to the true simulated DP assignments." << std::endl;
    std::cout << "   * -numSimulations    : The number of simulations to do inference under." << std::endl;
    std::cout << std::endl;

    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeMean          : Mean for the tree length prior." << std::endl;
    std::cout << "   * -treeSD            : SD for the tree length prior." << std::endl;
    std::cout << "   * -omegaLambda       : Rate parameter for the nonsynonymous mutation rate's exponential prior." << std::endl;
    std::cout << "   * -kLambda           : Rate parameter for the transition/transversion rate's exponential prior." << std::endl;
    std::cout << "   * -rLambda           : Rate parameter for the matrix-swapping rate's exponential prior." << std::endl;
    std::cout << "   * -expectedCat       : The number of expected categories for the DP." << std::endl;
    std::cout << "   * -fixedRegimes      : Fix the number of hidden regimes for CMM or DPCMM and disable RJ-MCMC over regime count." << std::endl;
    std::cout << "   * -mcmcmcBeta        : Inverse temperature for the auxiliary heated chain; values below 1.0 enable MCMCMC." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : The number of iterations for the MCMC." << std::endl;
    std::cout << "   * -printFreq         : How often to output the MCMC state to the screen." << std::endl;
    std::cout << "   * -sampleFreq        : How often to ouput the MCMC state to log files." << std::endl;
    std::cout << "   * -burnInIter        : The number of iterations for the burn-in." << std::endl;
    std::cout << "   * -tuneFreq          : How often to tune the acceptance rate of the MCMC moves during the burn-in." << std::endl;
    std::cout << "   * -bayesOpt          : (Experimental) The number of iterations to run Bayesian optimization on the MCMC moves after the burn-in." << std::endl;
    std::cout << "   * -bayesOptFreq      : (Experimental) How often to sample the trace for Bayesian optimization." << std::endl;
    std::cout << "   * -numGibbs          : How many Gibbs updates to perform on the DP clusters." << std::endl;
    std::cout << "   * -mcmcmcSwapFreq    : How many MCMC iterations between attempted swaps in MCMCMC." << std::endl;
    std::cout << "   * -treeWeight        : How often to propose a move on the tree." << std::endl;
    std::cout << "   * -kWeight           : How often to propose a move on the K parameter." << std::endl;
    std::cout << "   * -rWeight           : How often to propose a move on the R parameter." << std::endl;
    std::cout << "   * -rjWeight          : How often to propose a move on the dimension of the Markov modulated model." << std::endl;
    std::cout << "   * -stationaryWeight  : How often to propose a move on the stationary distribution." << std::endl;
    std::cout << "   * -dppWeight         : How often to propose a move on the DP clusters." << std::endl;
    std::cout << "   * -omegaWeight       : How often to propose a move on the omega parameters." << std::endl;
    std::cout << std::endl;
}
