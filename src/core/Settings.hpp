#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#include <string>
#include <vector>

struct Settings {
                    Settings(void) = delete;
                    Settings(int argc,  char* argv[]);
    std::string     nexusInput;
    std::string     treeOutput;
    std::string     dppOutput;
    std::string     mcmcOutput;
    std::string     tipsOutput;
    std::string     branchOutput;
    std::string     ancestralStatesOutput;
    std::string     simulationOutput;

    std::string     fixedTree;
    int             numIterations;
    int             printFrequency;
    int             sampleFrequency;
    int             burnInIterations;
    int             tuneFrequency;
    int             bayesOpt;
    int             bayesOptFrequency;
    int             truncation;
    double          kLambda;
    double          rLambda;
    double          omegaLambda;
    double          gammaLambda;
    double          expectedCat;
    double          treeLengthLambda;
    double          kWeight;
    double          rWeight;
    double          omegaWeight;
    double          dppWeight;
    double          treeWeight;
    double          proportionsWeight;
    double          stationaryWeight;
    double          rjWeight;
    int             numGibbs;
    int             numSimulations;

    bool M0;
    bool M3S2;
    bool SB;
    bool RJ;
    bool simulateM0;
    bool simulateM3S2;
    bool simulateDPP;
    bool simulateSB;
    bool sequentialTuningSim;

    void            usage();
    void            print();
};

#endif