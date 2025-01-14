#ifndef DIRICHLET_PROCESS_PRIOR_HPP
#define DIRICHLET_PROCESS_PRIOR_HPP
#include "modeling/parameters/Parameter.hpp"
#include <taskflow/taskflow.hpp>
#include <map>

class Model;
class Settings;

struct Category {
    double omega1;
    double omega2;
    double r;
    int size;
    std::vector<int> members;
    bool dirty;
};

class DirichletProcessPrior : public Parameter {
    public:
        DirichletProcessPrior(void)=delete;
        DirichletProcessPrior(int size, Settings s);
        ~DirichletProcessPrior(void);

        void registerModel(Model* m);
        
        double lnPrior() {return currentLnPrior;}

        double updateOmega();
        double updateR();
        double updateDPP();
        
        void tune();

        void accept();
        void reject();
        
        double getAlpha() {return alpha;}

        int getCategorySize(int index);
        void addCategory(double omega1, double omega2, double r);

        std::vector<Category> getCategories() {return currentCategories;}
        std::vector<int> getAssignments() { return assignments;}

        int unassignMember(int member);

        void assignMember(int member, int category);

        int getNumCategories(){return currentCategories.size();}

        int omegaAcceptCount;
        int omegaCount;
        int rAcceptCount;
        int rCount;
    protected:
        void regeneratePrior();
        Model* model;
        double currentLnPrior;
        double oldLnPrior;
        double alpha;
        double omegaLambda;
        double rLambda;

        double omegaDelta;
        double rDelta;

        int numMembers;
        int numGibbsUpdate;
        std::vector<Category> currentCategories;
        std::vector<Category> oldCategories;
        std::vector<int> assignments;

        int moveChoice;

        void removeCategory(int index);

        tf::Executor executor;
};

#endif