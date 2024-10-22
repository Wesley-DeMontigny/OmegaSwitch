#ifndef TRANSITION_PROBABILITY_HPP
#define TRANSITION_PROBABILITY_HPP

#include <complex>
#include <vector>
#include "core/EigenSystem.hpp"
#include "core/Matrix.hpp"
#include "RateEigens.hpp"

class TransitionProbability {

	public:
                                TransitionProbability(const int nn, const int nC);
                               ~TransitionProbability ();                                                                                     //!< destructor
        Matrix<double>*         operator()(int s, int r, int n) { return probs[s][r*numNodes + n]; }
        int                     getNumStates(void) { return numStates; }
        int                     getNumMatrices(void) {return isComplex.size();}
        void                    accept(void);          
        void                    reject(void);                                                                                   
        void                    setProbs(const int state, const int r, const int node, const double v);                                     //!< calculate transition probabilities (P) for length v and store it
        void                    updateQ(Matrix<double>& Q, const int index);
        void                    deleteQ(const int index);
        void                    popQ();

    private:
        EigenSystem*            eigens;
        std::vector<RateEigen>  rateEigen;
        std::vector<ComplexRateEigen> complexRateEigen;
        std::complex<double>*   ceigValExp;
        double*                 eigValExp;
        std::vector<bool>       isComplex;
        std::vector<bool>       isOldComplex;
        int                     numCats;
        int                     numNodes;
        int                     numStates;
        Matrix<double>**        probs[2];
        void                    calcCijk(int mIndex);
        void                    calcComplexCijk(int mIndex);
        void                    tiProbsComplexEigens(const double v, Matrix<double> &P, const int mIndex);
        void                    tiProbsEigens(const double v, Matrix<double> &P, const int mIndex);
};

#endif
