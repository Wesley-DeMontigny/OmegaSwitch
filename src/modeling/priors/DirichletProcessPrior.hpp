#ifndef DIRICHLET_PROCESS_PRIOR_HPP
#define DIRICHLET_PROCESS_PRIOR_HPP
#include "PriorNode.hpp"
#include "modeling/parameters/BasicParameter.hpp"

struct Category {
    double omega1;
    double omega2;
    int size;
    std::vector<int> members;
};

class DirichletProcessPrior : public PriorNode {
    public:
        DirichletProcessPrior(void)=delete;
        DirichletProcessPrior(int size, double alpha);
        ~DirichletProcessPrior(void);
        
        double lnPrior() {return currentLnPrior;}
        void regenerate();
        void accept();
        void reject();
        void sample();

        std::string writeValue() {return std::to_string(currentLnPrior);}
        
        double getAlpha() {return alpha;}

        int getCategorySize(int index);
        void addCategory(double value1, double value2);
        void setCategoryOmega1(int index, double value);
        void setCategoryOmega2(int index, double value);
        double getCategoryOmega1(int index);
        double getCategoryOmega2(int index);
        std::vector<Category> getCategories() {return currentCategories;}

        std::vector<int> getAssinments() { return assignments;}
        int unassignMember(int member);
        int unassignMember(int member, int category);
        int popBackCategory(int category);
        void assignMember(int member, int category);


        int getNumCategories(){return currentCategories.size();}
        int getNumMembers() {return numMembers;}
    protected:
        double oldLnPrior;
        double currentLnPrior;
        double alpha;
        int numMembers;
        std::vector<Category> currentCategories;
        std::vector<Category> oldCategories;
        std::vector<int> assignments;
        void removeCategory(int index);
};

#endif