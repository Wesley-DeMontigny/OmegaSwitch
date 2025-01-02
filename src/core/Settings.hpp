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
    int             numIterations;
    int             printFrequency;
    int             sampleFrequency;
    int             burnInIterations;
    int             tuneFrequency;
    double          rLambda;
    double          kLambda;
    double          omegaAlpha;
    double          omegaBeta;
    double          dppAlpha;
    double          treeLengthLambda;
    bool            updateStationary;
    int             numGibbsUpdate;
    double          rateMatrixWeight;
    double          dppWeight;
    double          treeWeight;

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
