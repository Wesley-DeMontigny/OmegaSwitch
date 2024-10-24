#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Settings.hpp"
#include "modeling/analysis/MoveScheduler.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/CodonMultiMatrix.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/analysis/Mcmc.hpp"
#include <chrono>

int main(int argc, char* argv[]) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    MoveScheduler moveScheduler;

    Settings settings(argc, argv);

    Alignment aln(settings.nexusInput);
    std::cout << "Initializing model..." << std::endl;

    std::vector<double> stationaryDist;
    for(double v : aln.getStateFrequencies()){
        stationaryDist.push_back(v/2);
    }
    for(double v : aln.getStateFrequencies()){
        stationaryDist.push_back(v/2);
    }

    TreeParameter treeParam(&aln, settings.treeLengthLambda);
    moveScheduler.registerParam(&treeParam, settings.treeWeight);

    DirichletProcessPrior dpp(aln.getNumChar(), settings.dppAlpha, settings.omegaLambda, settings.numGibbsUpdate);
    moveScheduler.registerParam(&dpp, settings.dppWeight);

    CodonMultiMatrix rateMatrix(settings.rLambda, settings.kLambda, stationaryDist, settings.updateStationary);
    moveScheduler.registerParam(&rateMatrix, settings.rateMatrixWeight);

    Model model(&aln, &treeParam, &rateMatrix, &dpp);

    Mcmc myMCMC(&model, &moveScheduler, settings);

    myMCMC.burnin();
    myMCMC.run();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
}