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
    std::string     ancestralStatesOutput;

    std::string     fixedTree;
    int             numIterations;
    int             printFrequency;
    int             sampleFrequency;
    int             burnInIterations;
    int             tuneFrequency;
    double          kLambda;
    double          rLambda;
    double          omegaLambda;
    double          expectedCat;
    double          treeLengthLambda;
    double          kWeight;
    double          rWeight;
    double          omegaWeight;
    double          dppWeight;
    double          treeWeight;
    double          stationaryWeight;
    int             numGibbs;

    void            usage();
    void            print();
};

#endif