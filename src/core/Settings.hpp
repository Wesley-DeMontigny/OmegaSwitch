#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#include <string>
#include <vector>

struct Settings {
    
                            Settings(int argc,  char* argv[]);
                            Settings(void) = delete;

    void                    print();
    void                    usage();

    bool                    M0 = false;
    bool                    RJ = false;
    bool                    RJDPP = false;
    bool                    sequentialTuningSim = false;
    bool                    simulateM0 = false;
    bool                    simulateRJ = false;
    bool                    simulateRJDPP = false;
    double                  dppWeight = 2.0;
    double                  expectedCat = 2.0;
    double                  kLambda = 5.0;
    double                  kWeight = 1.0;
    double                  omegaLambda = 5.0;
    double                  omegaWeight = 1.0;
    double                  rjWeight = 1.0;
    double                  rLambda = 5.0;
    double                  rWeight = 1.0;
    double                  stationaryWeight = 2.0;
    double                  treeLengthLambda = 5.0;
    double                  treeWeight = 2.0;
    int                     bayesOpt = 0;
    int                     bayesOptFrequency = 0;
    int                     burnInIterations = 1000;
    int                     numGibbs = 25;
    int                     numIterations = 5000;
    int                     numSimulations = 0;
    int                     printFrequency = 100;
    int                     sampleFrequency = 10;
    int                     tuneFrequency = 100;

    std::string             ancestralStatesOutput = "";
    std::string             branchOutput = "";
    std::string             dppOutput = "";
    std::string             fixedTree = "";
    std::string             mcmcOutput = "";
    std::string             nexusInput = "";
    std::string             simulationOutput = "";
    std::string             tipsOutput = "";
    std::string             treeOutput = "";
};

#endif