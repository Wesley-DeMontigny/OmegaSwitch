#ifndef RJ_MCMC_HPP
#define RJ_MCMC_HPP
#include <vector>
#include <string>
#include "core/BayesianOptimizer.hpp"

class Model;
class Move;
class Settings;

/**
 * @brief 
 * 
 */
class MCMC{
    public:
                                MCMC(void)=delete;                                              //
                                MCMC(Model* m, std::vector<Move>& mv, Settings& s, bool dBO);   //
        void                    burnin();                                                       //
        void                    run();                                                          //
    private:
        double                  GibbsIteration(double currentLnPosterior);                      //

        BayesianOptimizer       optim;                                                          //
        bool                    disableBayesOpt;                                                //
        double                  totalWeight;                                                    //
        int                     bayesOptFreq;                                                   //
        int                     bayesOptIter;                                                   //
        int                     numBurnIn;                                                      //
        int                     numIter;                                                        //
        int                     printFreq;                                                      //
        int                     sampleFreq;                                                     //
        int                     tuneFreq;                                                       //
        Model*                  model;                                                          //
        std::vector<Move>&      moves;                                                          //
};

#endif