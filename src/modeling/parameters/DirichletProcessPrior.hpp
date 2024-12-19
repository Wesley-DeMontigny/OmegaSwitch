#ifndef DIRICHLET_PROCESS_PRIOR_HPP
#define DIRICHLET_PROCESS_PRIOR_HPP
#include "modeling/parameters/Parameter.hpp"
#include <taskflow/taskflow.hpp>
#include <map>

class Model;

struct Category {
    double value;
    int size;
    std::vector<int> members;
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
        void addCategory(double omega);
        std::vector<int> getAssinments() { return assignments;}

        int unassignMember(int member);

        void assignMember(int member, int category);

        int getNumCategories(){return categories.size();}
    protected:
        void regeneratePrior();
        Model* model;
        double currentLnPrior;
        double alpha;
        double omegaLambda;

        int numMembers;
        int numGibbsUpdate;
        std::vector<Category> categories;
        std::vector<double> assignments;

        void removeCategory(int index);

        tf::Executor executor;
};

#endif