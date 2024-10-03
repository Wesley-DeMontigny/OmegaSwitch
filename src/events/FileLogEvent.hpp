#ifndef FILE_LOG_EVENT
#define FILE_LOG_EVENT
#include "Event.hpp"
#include <string>
#include <vector>

class ModelNode;

class FileLogEvent : public Event{
    public:
        FileLogEvent(void)=delete;
        FileLogEvent(std::vector<std::pair<std::string, ModelNode*>> n, std::string f);
        void initialize();
        void call(int iteration);
    private:
        std::vector<std::pair<std::string, ModelNode*>> nodes;
        std::string file;
};

#endif