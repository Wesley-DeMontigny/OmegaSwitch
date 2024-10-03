#include "FileLogEvent.hpp"
#include "modeling/ModelNode.hpp"
#include "core/Msg.hpp"
#include <string>
#include <vector>
#include <fstream>

FileLogEvent::FileLogEvent(std::vector<std::pair<std::string, ModelNode*>> n, std::string f) : nodes(n), file(f) {}

void FileLogEvent::initialize() {
    std::fstream fs;
    fs.open (file, std::fstream::out);

    fs << "Iteration";

    for(std::pair<std::string, ModelNode*> entry : nodes)
        fs << "\t" << entry.first;
    fs << "\n";

    fs.close();
}

void FileLogEvent::call(int iteration) {
    std::fstream fs;
    fs.open (file, std::fstream::app);

    fs << iteration;

    for(std::pair<std::string, ModelNode*> entry : nodes)
        fs << "\t" << entry.second->writeValue();
    fs << "\n";

    fs.close();
}