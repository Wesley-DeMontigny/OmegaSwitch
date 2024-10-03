#include "IterationTrackerEvent.hpp"
#include <iostream>

IterationTrackerEvent::IterationTrackerEvent(void){}

void IterationTrackerEvent::initialize() {}

void IterationTrackerEvent::call(int iteration) {
    std::cout << "Iteration " << iteration << " Completed..." << std::endl;
}