#include "DPPFileLogEvent.hpp"
#include "modeling/PosteriorNode.hpp"
#include "core/Msg.hpp"
#include <string>
#include <vector>
#include <fstream>

DPPFileLogEvent::DPPFileLogEvent(DirichletProcessPrior* d, PosteriorNode* p, std::string f, bool normalize) : dpp(d), posterior(p), file(f), normalizeCategoryValues(normalize) {}

void DPPFileLogEvent::initialize() {
    std::fstream fs;
    fs.open (file, std::fstream::out);

    fs << "Iteration\tPosterior\tCategories";

    for(int i = 0, len = dpp->getNumMembers(); i < len; i++)
        fs << "\tDPP[" << i << "]";
    fs << "\n";

    fs.close();
}

void DPPFileLogEvent::call(int iteration) {
    std::fstream fs;
    fs.open (file, std::fstream::app);

    fs << iteration << "\t" << posterior->lnPosterior() << "\t" << dpp->getNumCategories();

    std::vector<int> assignments = dpp->getAssinments();
    std::vector<double> catValues = dpp->getCategoryValues();
    if(normalizeCategoryValues){
        double maxL = *std::max_element(catValues.begin(), catValues.end());
        for(double& cat : catValues)
            cat = cat/maxL;
    }

    for(int a : assignments)
        fs << "\t" << catValues[a];
    fs << "\n";

    fs.close();
}