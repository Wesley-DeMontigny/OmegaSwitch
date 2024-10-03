#ifndef DPP_FILE_LOG_EVENT
#define DPP_FILE_LOG_EVENT
#include "Event.hpp"
#include <string>
#include <vector>
#include "modeling/priors/DirichletProcessPrior.hpp"

class PosteriorNode;

class DPPFileLogEvent : public Event{
    public:
        DPPFileLogEvent(void)=delete;
        DPPFileLogEvent(DirichletProcessPrior* d, PosteriorNode* p, std::string f);
        void initialize();
        void call(int iteration);
    private:
        DirichletProcessPrior* dpp;
        PosteriorNode* posterior;
        std::string file;
        bool normalizeCategoryValues;
};

#endif