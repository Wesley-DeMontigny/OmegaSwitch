#include "ScreenLogEvent.hpp"
#include "modeling/ModelNode.hpp"
#include "core/Msg.hpp"
#include <string>
#include <vector>
#include <iostream>

ScreenLogEvent::ScreenLogEvent(std::vector<std::pair<std::string, ModelNode*>> n) : nodes(n) {}

void ScreenLogEvent::initialize() {
    std::cout << "Iteration";

    for(std::pair<std::string, ModelNode*> entry : nodes)
        std::cout << "\t" << entry.first;
    std::cout << std::endl;
}

void ScreenLogEvent::call(int iteration) {
    std::cout << iteration;

    for(std::pair<std::string, ModelNode*> entry : nodes)
        std::cout << "\t" << entry.second->writeValue();
    std::cout << std::endl;
}