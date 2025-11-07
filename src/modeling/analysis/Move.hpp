#ifndef MOVE_HPP
#define MOVE_HPP
#include <functional>
#include <vector>
#include <optional>
#include <iostream>

/**
 * @brief An MCMC "move" or proposal that is called to update the parameters
 * in a particular way.
 */
struct Move {
    double                      weight;                         // The weight of the proposal (how often to propose it)
    int                         gibbsIterations;                // The number of Metropolis-Hastings proposals per Gibbs sampling iteration
    std::function<double()>     action;                         // The funciton that calls the move
    std::function<bool()>       condition;                      // The condition for including this move
};

#endif