#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#include <string>

struct Settings {
                    Settings(void) = delete;
                    Settings(int argc,  char* argv[]);
    std::string     nexusInput;
    std::string     treeOutput;
    std::string     dppOutput;
    std::string     mcmcOutput;
    int             numIterations;
    int             printFrequency;
    int             sampleFrequency;
    int             burnInIterations;
    int             tuneFrequency;
    double          rLambda;
    double          kLambda;
    double          omegaLambda;
    double          dppAlpha;
    double          treeLengthLambda;
    bool            updateStationary;
    int             numGibbsUpdate;
    double          rateMatrixWeight;
    double          dppWeight;
    double          treeWeight;

    void usage();
    void print();
};

#endif
