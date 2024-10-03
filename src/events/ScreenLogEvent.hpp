#ifndef SCREEN_LOG_EVENT
#define SCREEN_LOG_EVENT
#include "Event.hpp"
#include <string>
#include <vector>

class ModelNode;

class ScreenLogEvent : public Event{
    public:
        ScreenLogEvent(void)=delete;
        ScreenLogEvent(std::vector<std::pair<std::string, ModelNode*>> n);
        void initialize();
        void call(int iteration);
    private:
        std::vector<std::pair<std::string, ModelNode*>> nodes;
};

#endif