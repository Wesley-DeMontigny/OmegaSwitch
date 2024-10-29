#ifndef DIRICHLET_PROCESS_PRIOR_HPP
#define DIRICHLET_PROCESS_PRIOR_HPP
#include "modeling/parameters/Parameter.hpp"
#include <map>

class Model;

struct Category {
    double omega;
    double beta;
    int size;
    std::vector<int> members;
    bool dirty;
};

class DirichletProcessPrior : public Parameter {
    public:
        DirichletProcessPrior(void)=delete;
        DirichletProcessPrior(int size, double alpha, double oL, int numGibbs);
        ~DirichletProcessPrior(void);

        void registerModel(Model* m) {model = m;}
        
        double lnPrior() {return currentLnPrior;}
        double update();
        void tune();

        void accept();
        void reject();
        
        double getAlpha() {return alpha;}

        int getCategorySize(int index);
        void addCategory(double value1, double value2);
        double getCategoryOmega(int index) {return currentCategories[index].omega;}
        double getCategoryBeta(int index) {return currentCategories[index].beta;}
        std::vector<Category> getCategories() {return currentCategories;}
        Category getCategory(int index) {return currentCategories[index];}

        std::vector<int> getAssinments() { return assignments;}
        int unassignMember(int member);
        int unassignMember(int member, int category);
        int popBackCategory(int category);
        void assignMember(int member, int category);

        int getNumCategories(){return currentCategories.size();}
        int getNumMembers() {return numMembers;}
    protected:
        void regeneratePrior();
        Model* model;
        double oldLnPrior;
        double currentLnPrior;
        double alpha;
        double omegaLambda;

        int moveChoice;
        int omegaCount;
        int omegaAcceptCount;
        double omegaDelta;
        int betaCount;
        int betaAcceptCount;
        double betaAlpha;

        int numMembers;
        int numGibbsUpdate;
        double denominator;
        std::vector<Category> currentCategories;
        std::vector<Category> oldCategories;
        std::vector<int> assignments;
        std::vector<int> oldAssignments;

        void removeCategory(int index);
};

#endif