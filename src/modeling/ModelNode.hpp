#ifndef MODEL_NODE_HPP
#define MODEL_NODE_HPP
#include <vector>
#include <string>

class ModelNode {
    public:
        virtual void accept()=0;
        virtual void reject()=0;
        virtual void regenerate()=0;
        virtual std::string writeValue()=0;
        void clean() {dirtyFlag = false;}
        void dirty() {dirtyFlag = true;}
        bool isDirty() {return dirtyFlag;}
    private:
        bool dirtyFlag;
};

#endif