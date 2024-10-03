#ifndef RATE_MATRIX_HPP
#define RATE_MATRIX_HPP
#include "modeling/ModelNode.hpp"
#include "core/Matrix.hpp"
#include <vector>

class PriorNode;
class LikelihoodNode;

class RateMatrix : public ModelNode {
    public:
        virtual Matrix<double> Q()=0;
        virtual std::vector<double> stationary()=0;
        virtual void accept()=0;
        virtual void reject()=0;
        virtual void regenerate()=0;
        virtual std::string writeValue()=0;
};

#endif