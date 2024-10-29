#ifndef PARAMETER_HPP
#define PARAMETER_HPP
#include <vector>
#include <string>

class Parameter {
    public:
        virtual void accept()=0;
        virtual void reject()=0;
        virtual double lnPrior()=0;
        virtual double update()=0;
        virtual void tune()=0;
        void clean() {dirtyFlag = false;}
        void dirty() {dirtyFlag = true;}
        bool isDirty() {return dirtyFlag;}
    private:
        bool dirtyFlag;
};

#endif