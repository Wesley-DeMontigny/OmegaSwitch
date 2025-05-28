#ifndef TRANSITION_PROBABILITY_HPP
#define TRANSITION_PROBABILITY_HPP

#include <complex>
#include <vector>
#include "core/EigenSystem.hpp"
#include "core/Matrix.hpp"
#include "core/RateEigens.hpp"

class TransitionProbability {

	public:
                                TransitionProbability(const int nn, const int ss);
                               ~TransitionProbability ();
        const Matrix<double>&   operator()(int s, int r, int n) const {
                                    return (s == 0) ? probs1[r][n] : probs2[r][n];
                                }
        int                     getNumStates(void) { return numStates; }
        int                     getNumMatrices(void) {return isComplex.size();}
        void                    accept(void);          
        void                    reject(void);                                                                                   
        void                    setProbs(const int state, const int r, const int node, const double v);
        std::vector<Matrix<double>>     generateProbs(Matrix<double> Q, std::vector<double> branches);
        void                    updateQ(Matrix<double> Q, const int index);
        void                    deleteQ(const int index);
        void                    deleteNQ(const int count);
        void                    allocateQ(int size);

    private:
        EigenSystem*            eigens;
        std::vector<RateEigen>  rateEigen;
        std::vector<ComplexRateEigen> complexRateEigen;
        std::vector<bool>       isComplex;
        std::vector<bool>       isOldComplex;
        int                     numNodes;
        int                     numStates;
        std::vector<Matrix<double>*> probs1;
        std::vector<Matrix<double>*> probs2;
        void                    tiProbsComplexEigens(const double v, Matrix<double> &P, ComplexRateEigen& rE);
        void                    tiProbsEigens(const double v, Matrix<double> &P, RateEigen& rE);
};

#endif
