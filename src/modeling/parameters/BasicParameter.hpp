#ifndef BASIC_PARAMETER_HPP
#define BASIC_PARAMETER_HPP
#include "modeling/ModelNode.hpp"

template <class T>
class BasicParameter : public ModelNode {
    public:
        BasicParameter(void)=delete;
        BasicParameter(T v) : currentValue(v), oldValue(v) {}
        void accept() {oldValue = currentValue;}
        void reject() {currentValue = oldValue;}
        void regenerate() {};
        double getValue() {return currentValue;}
        void setValue(T v) {currentValue = v;}
        std::string writeValue() {return std::to_string(currentValue);}
    protected:
        T currentValue;
        T oldValue;
};

#endif