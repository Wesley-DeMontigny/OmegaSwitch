#include "moves/MoveScheduler.hpp"
#include "TuneEvent.hpp"
#include "core/Msg.hpp"

TuneEvent::TuneEvent(MoveScheduler* m) : moveScheduler(m) {}

void TuneEvent::initialize() {
    if(moveScheduler == nullptr)
        Msg::error("TuneEvent expects moveScheduler to be set!");
}

void TuneEvent::call(int iteration) {
    moveScheduler->tune();
}