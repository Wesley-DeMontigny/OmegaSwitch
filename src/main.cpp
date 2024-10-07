#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "moves/MoveScheduler.hpp"
#include "moves/Move.hpp"
#include "moves/MoveTreeLocal.hpp"
#include "moves/MoveScaleBranch.hpp"
#include "moves/MoveTreeNNI.hpp"
#include "moves/MoveBetaSimplex.hpp"
#include "moves/MoveDPPCodonGibbs.hpp"
#include "moves/MoveScaleDouble.hpp"
#include "moves/MoveScaleDPPCategory.hpp"
#include "events/EventManager.hpp"
#include "events/TuneEvent.hpp"
#include "events/FileLogEvent.hpp"
#include "events/ScreenLogEvent.hpp"
#include "events/IterationTrackerEvent.hpp"
#include "events/DPPFileLogEvent.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/rates/CodonMultiMatrix.hpp"
#include "modeling/parameters/BasicParameter.hpp"
#include "modeling/likelihoods/PhyloCTMC.hpp"
#include "modeling/priors/TreePrior.hpp"
#include "modeling/priors/DirichletPrior.hpp"
#include "modeling/priors/DirichletProcessPrior.hpp"
#include "modeling/priors/ExponentialRatioPrior.hpp"
#include "modeling/PosteriorNode.hpp"
#include "modeling/analysis/Mcmc.hpp"
#include "modeling/priors/GammaPrior.hpp"
#include <chrono>

int main(int argc, char* argv[]) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    MoveScheduler moveScheduler;

    Alignment aln("C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/replicase.nex");

    //Purely empirical values right now.
    std::vector<BasicParameter<double>*> stationaryDist;
    for(double v : aln.getStateFrequencies()){
        stationaryDist.push_back(new BasicParameter<double>(v/2));
    }
    for(double v : aln.getStateFrequencies()){
        stationaryDist.push_back(new BasicParameter<double>(v/2));
    }

    TreeParameter treeParam(&aln);
    TreePrior treePrior(&treeParam);
    BasicParameter<double> lambda(10.0);
    treePrior.setExponentialBranchPrior(&lambda);
    treePrior.sample();

    BasicParameter<double> k(1.0);
    ExponentialRatioPrior kPrior(&k);
    kPrior.sample();
    MoveScaleDouble kMove(&k);
    moveScheduler.registerMove(&kMove, 2.5);

    BasicParameter<double> r(1.0);
    ExponentialRatioPrior rPrior(&r);
    rPrior.sample();
    MoveScaleDouble rMove(&r);
    moveScheduler.registerMove(&rMove, 2.5);

    DirichletProcessPrior dpp(aln.getNumChar(), 0.5);

    CodonMultiMatrix rateMatrix(&dpp, &k, &r, stationaryDist);

    PhyloCTMC ctmc(&aln, &treeParam, &rateMatrix, &dpp);

    MoveDPPCodonGibbs moveDPP(&ctmc, &dpp);
    moveScheduler.registerMove(&moveDPP, 10.0);
    MoveScaleDPPCategory moveCategories(&dpp);
    moveScheduler.registerMove(&moveCategories, 10.0);

    MoveTreeLocal localMove(&treeParam);
    moveScheduler.registerMove(&localMove, 15.0);

    PosteriorNode posterior(&ctmc, {&treePrior, &kPrior, &rPrior, &dpp});

    Mcmc myMCMC(&posterior, &moveScheduler);

    EventManager burnIn;
    burnIn.registerEvent(&TuneEvent(&moveScheduler), 250);
    burnIn.registerEvent(&IterationTrackerEvent(), 250);
    burnIn.initialize();

    std::cout << "Starting Burn-In" << std::endl;
    //myMCMC.run(30000, &burnIn);

    std::vector<std::pair<std::string, ModelNode*>> loggables;
    loggables.push_back(std::make_pair("Tree Prior", &treePrior));
    loggables.push_back(std::make_pair("K Prior", &kPrior));
    loggables.push_back(std::make_pair("R Prior", &rPrior));
    loggables.push_back(std::make_pair("DPP Prior", &dpp));
    loggables.push_back(std::make_pair("Likelihood", &ctmc));
    loggables.push_back(std::make_pair("Posterior", &posterior));
    loggables.push_back(std::make_pair("K", &k));
    loggables.push_back(std::make_pair("R", &r));

    EventManager realRun;
    ScreenLogEvent screenLogger(loggables);
    realRun.registerEvent(&screenLogger, 1);
    FileLogEvent fileLogger(loggables, "C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/test_mcmc.log");
    realRun.registerEvent(&fileLogger, 1);
    loggables.push_back(std::make_pair("Tree", &treeParam));
    loggables.erase(loggables.begin());
    loggables.erase(loggables.begin());
    FileLogEvent treeLogger(loggables, "C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/tree_trace.trees");
    realRun.registerEvent(&treeLogger, 1);

    DPPFileLogEvent logDPP(&dpp, &posterior, "C:/Users/wescd/OneDrive/Documents/Code/Varying_Selection_DPP/res/dpp_sites.log");
    realRun.registerEvent(&logDPP, 1);

    realRun.initialize();

    myMCMC.run(250000, &realRun);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    std::cout << treeParam.getTree()->getNewick() << std::endl;

    for(BasicParameter<double>* f : stationaryDist)
        delete f;
}
