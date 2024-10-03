#include "DPPFileLogEvent.hpp"
#include "modeling/PosteriorNode.hpp"
#include "core/Msg.hpp"
#include <string>
#include <vector>
#include <fstream>

DPPFileLogEvent::DPPFileLogEvent(DirichletProcessPrior* d, PosteriorNode* p, std::string f) : dpp(d), posterior(p), file(f) {}

void DPPFileLogEvent::initialize() {
    std::fstream fs;
    fs.open (file, std::fstream::out);

    fs << "Iteration\tPosterior\tCategories";

    for(int i = 0, len = dpp->getNumMembers(); i < len; i++)
        fs << "\tDPP[" << i << ",1]" << "\tDPP[" << i << ",2]";
    fs << "\n";

    fs.close();
}

void DPPFileLogEvent::call(int iteration) {
    std::fstream fs;
    fs.open (file, std::fstream::app);

    fs << iteration << "\t" << posterior->lnPosterior() << "\t" << dpp->getNumCategories();

    std::vector<int> assignments = dpp->getAssinments();
    std::vector<Category> categories = dpp->getCategories();

    for(int a : assignments)
        fs << "\t" << categories[a].omega1 << "\t" << categories[a].omega2;
    fs << "\n";

    fs.close();
}