#ifndef TRANSITION_PROBABILITY_HPP
#define TRANSITION_PROBABILITY_HPP

#include <complex>
#include <vector>
#include "core/EigenSystem.hpp"
#include "core/Matrix.hpp"
#include "RateEigens.hpp"

/* This class sets up a continuous-time Markov model describing how discrete characters (DNA, RNA, amino acid,
   morphological etc.) change on a phylogenetic tree. A continuous-time Markov chain is defined by a matrix describing
   the infinitessimal rates of change from one state to another. This rate matrix, also known as the Q matrix or the
   instantaneous rate matrix, allows one to calculate the probability of observing a change from state i to state
   j over a branch length (evolutionary time) of v. The probabilities of change can be expressed in a matrix known
   as the P matrix or the transition probability matrix. The instantaneous rate matrix also allows one to determine
   the long-term behavior of the chain, that is, the probability of finding the process in state i after infinite
   time. This is also known as the stationary frequency of i.
  
   The class accommodates both reversible and irreversible models. An irreversible model is created by a constructor
   taking a rate set parameter and a data type parameter (defaulting to DNA). A reversible model is created by a
   separate constructor, which also takes a state frequency parameter. Both constructors also take a bool indicating
   whether the Eigensystem or the Pade approximation will be used to calculate the Q matrix from the P matrix.
  
   The rates in an reversible model are assumed to be in the following order (illustrated for the DNA data type):
  
            to
             A      C      G      T
   from A  -----  r0*f1  r1*f2  r2*f3
        C  r0*f0  -----  r3*f2  r4*f3
        G  r1*f0  r3*f1  -----  r5*f3
        T  r2*f0  r4*f1  r5*f2  -----
  
   That is, the rates are in the order of the upper diagonal: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T.
  
   For an irreversible model, the rates are assumed to be in the order (for a DNA character):
  
            to
            A    C    G    T
   from A  ---  r0   r1   r2
        C  r3   ---  r4   r5
        G  r6   r7   ---  r8
        T  r9   r10  r11  ---
  
   That is, the off-diagonal rates are given one row at a time: A->C, A->G, A->T, C->A, C->G,
   C->T, G->A, G->C, G->T, T->A, T->C, T->G. */

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
        int                     updateQ(const Matrix<double>& qTemp, const int index);
        void                    deleteQ(const int index);
        void                    popQ();

    private:
        EigenSystem*            eigens;
        std::vector<RateEigen>  rateEigen;
        std::vector<ComplexRateEigen> complexRateEigen;
        std::complex<double>*   ceigValExp;                                                                                                //!< working space for calculating exp(-lambda*v) from complex eigensystem
        double*                 eigValExp;                                                                                                 //!< working space for calculating exp(-lambda*v)
        std::vector<bool>       isComplex;                                                                                                  //!< does Q have complex eigensystem?
        std::vector<bool>       isOldComplex;
        int                     numCats;
        int                     numNodes;
        int                     numStates;                                                                                                  //!< number of states                                                                                      
        Matrix<double>          Q;                                                                                                          //!< the Q (rate) matrix
        std::vector<double>     pi;                                                                                                       //!< stationary frequencies (for irrev matrix)
        Matrix<double>**        probs[2];
        void                    calcCijk(int mIndex);                                                                                             //!< precalculations for matrix exponentiation using eigensystem
        void                    calcComplexCijk(int mIndex);                                                                                      //!< precalculations for matrix exponentiation using complex eigensystem
        void                    calcStationaryFreq(void);                                                                                   //!< calculate the stationary probabilites
        void                    initializeEigenVariables(void);                                                                             //!< initialize local variables for eigensystem calculations
        void                    initializeProbabilityBuffer();
        void                    rescaleQ(bool scalePi = false);                                                                                             //!< rescale Q matrix (using pi)
        void                    tiProbsComplexEigens(const double v, Matrix<double> &P, const int mIndex);                                                  //!< calculates transition probabilities using complex eigensystem
        void                    tiProbsEigens(const double v, Matrix<double> &P, const int mIndex);                                                         //!< calculates transition probabilities using eigensystem
};

#endif
