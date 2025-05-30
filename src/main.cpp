#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/parameters/DirichletProcessPrior.hpp"
#include "modeling/parameters/DPPMatrix.hpp"
#include "modeling/model/DPPModel.hpp"
#include "modeling/analysis/DPPMcmc.hpp"
#include "modeling/parameters/M0Matrix.hpp"
#include "modeling/model/M0Model.hpp"
#include "modeling/analysis/M0Mcmc.hpp"
#include "modeling/parameters/M3S2Matrix.hpp"
#include "modeling/model/M3S2Model.hpp"
#include "modeling/analysis/M3S2Mcmc.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <algorithm>
#include <chrono>

void inference(Settings settings, Alignment aln, TreeParameter treeParam){
    if(settings.M0){
        std::cout << "Initializing the M0 model..." << std::endl;

        M0Matrix rateMatrix(settings);

        M0Model model(settings, &aln, &treeParam, &rateMatrix);

        M0Mcmc myMCMC(&model, &treeParam, &rateMatrix, settings);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else if(settings.M3S2){
        std::cout << "Initializing the M3S2 model..." << std::endl;

        M3S2Matrix rateMatrix(settings);

        M3S2Model model(settings, &aln, &treeParam, &rateMatrix);

        M3S2Mcmc myMCMC(&model, &treeParam, &rateMatrix, settings);
        
        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
    else {
        std::cout << "Initializing the DPP model..." << std::endl;
        
        DirichletProcessPrior dpp(aln.getNumChar(), settings);

        DPPMatrix rateMatrix(settings);

        DPPModel model(settings, &aln, &treeParam, &rateMatrix, &dpp);

        DPPMcmc myMCMC(&model, &treeParam, &rateMatrix, &dpp, settings);

        std::cout << "Starting MCMC..." << std::endl;
        myMCMC.burnin();
        myMCMC.run();
    }
}

int main(int argc, char* argv[]) {

    Settings settings(argc, argv);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(settings.simulateDPP == false && settings.simulateM0 == false && settings.simulateM3S2 == false){
        Alignment aln(settings.nexusInput);
        TreeParameter treeParam(&aln, settings.fixedTree, settings.treeLengthLambda);
        inference(settings, aln, treeParam);
    }
    else if(settings.simulateDPP){
        // TODO
    }
    else if(settings.simulateM0){
        // TODO
    }
    else if(settings.simulateM3S2){
        // TODO
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
}