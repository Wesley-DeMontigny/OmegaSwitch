#ifndef RJ_DIRICHLET_PROCESS_PRIOR_HPP
#define RJ_DIRICHLET_PROCESS_PRIOR_HPP
#include "modeling/parameters/Parameter.hpp"
#include <taskflow/taskflow.hpp>
#include <map>

class RJDPPModel;
class Settings;
class RJDPPMatrix;

struct RJCategory {
    double omega1;
    double omega2;
    double omega3;
    double omegaIncrement1;
    double omegaIncrement2;
    int size;
    std::vector<int> members;
    bool dirty;
};

class RJDirichletProcessPrior : public Parameter {
    public:
        RJDirichletProcessPrior(void)=delete;
        RJDirichletProcessPrior(RJDPPMatrix* matrix, int size, Settings s);
        ~RJDirichletProcessPrior(void);

        void registerModel(RJDPPModel* m);
        
        double lnPrior() {return currentLnPrior;}

        double updateOmega();
        double updateDPP();
        double updateActiveOmegas();
        
        void tune();

        void accept();
        void reject();
        
        double getAlpha() {return alpha;}

        std::vector<RJCategory> getCategories() {return currentCategories;}
        std::vector<int> getAssignments() { return assignments;}
        double getActiveOmegas() {return currentActiveOmegas;}

        int getNumCategories(){return currentCategories.size();}

        double omegaDelta;
        int omegaAcceptCount;
        int omegaCount;
    private:
        double calculateAlpha(double expectedCat, int members);
        double expectedCategories(double a, int members);
        void regeneratePrior();

        int getCategorySize(int index);
        void addCategory(double omega1, double omega2, double omega3);

        int unassignMember(int member);
        void assignMember(int member, int category);

        RJDPPModel* model;
        RJDPPMatrix* rateMatrix;

        double currentLnPrior;
        double oldLnPrior;

        int currentActiveOmegas = 3;
        int oldActiveOmegas = 3;

        double alpha;
        double omegaLambda;

        int numMembers;
        int numGibbs;
        std::vector<RJCategory> currentCategories;
        std::vector<RJCategory> oldCategories;
        std::vector<int> assignments;

        int moveChoice;

        void removeCategory(int index);

        tf::Executor executor;
};

#endif