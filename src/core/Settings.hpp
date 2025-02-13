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

    std::string     fixedTree;
    bool            normalizeRates;
    int             numIterations;
    int             printFrequency;
    int             sampleFrequency;
    int             burnInIterations;
    int             tuneFrequency;
    double          kLambda;
    double          rLambda;
    double          omegaLambda;
    double          dppAlpha;
    double          treeLengthLambda;
    int             numGibbsUpdate;
    double          kWeight;
    double          rWeight;
    double          omegaWeight;
    double          dppWeight;
    double          treeWeight;
    double          stationaryWeight;

    bool            simulate;
    int             numTaxa;
    int             numChar;
    double          kValue;
    double          rValue;
    std::vector<int> assignmentVector;
    std::vector<double> omega1Vector;
    std::vector<double> omega2Vector;

    void            usage();
    void            print();
};

#endif
