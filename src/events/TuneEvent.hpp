#ifndef TUNE_EVENT_HPP
#define TUNE_EVENT_HPP
#include "Event.hpp"

class MoveScheduler;

class TuneEvent : public Event{
    public:
        TuneEvent(void)=delete;
        TuneEvent(MoveScheduler* m);
        void initialize();
        void call(int iteration);
    private:
        MoveScheduler* moveScheduler;
};

#endif;