#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/DPPMatrix.hpp"
#include "modeling/model/DPPModel.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/analysis/DPPMcmc.hpp"
#include "modeling/parameters/M0Matrix.hpp"
#include "modeling/model/M0Model.hpp"
#include "modeling/analysis/M0Mcmc.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <algorithm>
#include <chrono>

int main(int argc, char* argv[]) {

    Settings settings(argc, argv);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Alignment aln(settings.nexusInput);
    std::cout << "Initializing model..." << std::endl;

    TreeParameter treeParam(&aln, settings.fixedTree, settings.treeLengthLambda);
    
    DirichletProcessPrior dpp(aln.getNumChar(), settings);

    DPPMatrix rateMatrix(settings);

    DPPModel model(settings, &aln, &treeParam, &rateMatrix, &dpp);

    DPPMcmc myMCMC(&model, &treeParam, &rateMatrix, &dpp, settings);

    /*
    M0Matrix rateMatrix(settings);

    M0Model model(settings, &aln, &treeParam, &rateMatrix);

    M0Mcmc myMCMC(&model, &treeParam, &rateMatrix, settings);
    */

    std::cout << "Starting MCMC..." << std::endl;
    myMCMC.burnin();
    myMCMC.run();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
}