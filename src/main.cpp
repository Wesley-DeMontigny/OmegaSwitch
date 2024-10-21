#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
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

    Alignment aln("C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/replicase.nex");

    //Purely empirical values right now.
    std::vector<double> stationaryDist;
    for(double v : aln.getStateFrequencies()){
        stationaryDist.push_back(v/2);
    }
    for(double v : aln.getStateFrequencies()){
        stationaryDist.push_back(v/2);
    }

    TreeParameter treeParam(&aln, 10.0);
    moveScheduler.registerParam(&treeParam, 20.0);

    DirichletProcessPrior dpp(aln.getNumChar(), 0.5, 5);
    moveScheduler.registerParam(&dpp, 10.0);

    CodonMultiMatrix rateMatrix(2.0, stationaryDist, false);
    moveScheduler.registerParam(&rateMatrix, 10.0);

    Model model(&aln, &treeParam, &rateMatrix, &dpp);

    Mcmc myMCMC(&model, &moveScheduler);

    myMCMC.burnin(25000, 10, 5000);
    myMCMC.run(500000, 10, 100);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
}