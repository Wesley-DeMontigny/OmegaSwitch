#ifndef MOVE_HPP
#define MOVE_HPP
#include <functional>
#include <vector>
#include <optional>
#include <iostream>

/**
 * @brief 
 * 
 */
struct Move {
    double                      weight;                         //
    int                         gibbsIterations;                //
    std::function<double()>     action;                         //
    std::function<bool()>       condition;                      //
};

#endif