#ifndef ITERATION_TRACKER_EVENT
#define ITERATION_TRACKER_EVENT
#include "Event.hpp"

class LikelihoodNode;
class PriorNode;
class ParameterNode;

class IterationTrackerEvent : public Event{
    public:
        IterationTrackerEvent(void);
        void initialize();
        void call(int iteration);
    private:
};

#endif