#ifndef RATE_EIGEN_HPP
#define RATE_EIGEN_HPP
#include <complex>
#include <vector>



/**
 * @brief Struct used to contain all of the values needed to do rate matrix exponentiation
 * via eigen decomposition when the eigenvalues are all real. This also contains all information
 * needed to revert to the old values if there is an update and a rejection during MCMC
 */
struct RateEigen {
    double*     c_ijk;          // This is a pre-computed tensor to avoid computing V D V^-1 every time. Cijk = Vik (V^-1)kj
    double*     eigenvalue;     // The real eigenvalues obtained by EigenSystem
    double*     oldC_ijk;       // Old values to reset the current ones if there was a rejection
    double*     oldEigenvalue;
    int         numStates;

    // Constructor
    RateEigen(int nS) : numStates(nS) {
        c_ijk = new double[numStates * numStates * numStates];
        eigenvalue = new double[numStates];
        oldC_ijk = new double[numStates * numStates * numStates];
        oldEigenvalue = new double[numStates];

        std::fill(c_ijk, c_ijk + numStates * numStates * numStates, 0.0);
        std::fill(oldC_ijk, oldC_ijk + numStates * numStates * numStates, 0.0);
        std::fill(eigenvalue, eigenvalue + numStates, 0.0);
        std::fill(oldEigenvalue, oldEigenvalue + numStates, 0.0);
    }

    // Copy constructor
    RateEigen(const RateEigen& other) {
        numStates = other.numStates;
        c_ijk = new double[numStates * numStates * numStates];
        eigenvalue = new double[numStates];
        oldC_ijk = new double[numStates * numStates * numStates];
        oldEigenvalue = new double[numStates];
        
        std::copy(other.c_ijk, other.c_ijk + numStates * numStates * numStates, c_ijk);
        std::copy(other.oldC_ijk, other.oldC_ijk + numStates * numStates * numStates, oldC_ijk);
        std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
        std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
    }

    // Assignment operator
    RateEigen& operator=(const RateEigen& other) {
        if (this != &other) {
            delete [] c_ijk;
            delete [] eigenvalue;
            delete [] oldC_ijk;
            delete [] oldEigenvalue;

            numStates = other.numStates;
            c_ijk = new double[numStates * numStates * numStates];
            eigenvalue = new double[numStates];
            oldC_ijk = new double[numStates * numStates * numStates];
            oldEigenvalue = new double[numStates];

            std::copy(other.c_ijk, other.c_ijk + numStates * numStates * numStates, c_ijk);
            std::copy(other.oldC_ijk, other.oldC_ijk + numStates * numStates * numStates, oldC_ijk);
            std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
            std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
        }
        return *this;
    }

    ~RateEigen() {
        delete [] c_ijk;
        delete [] eigenvalue;
        delete [] oldEigenvalue;
        delete [] oldC_ijk;
    }
};

/**
 * @brief Struct used to contain all of the values needed to do rate matrix exponentiation
 * via eigen decomposition when the eigenvalues are complex. This also contains all information
 * needed to revert to the old values if there is an update and a rejection during MCMC
 */
struct ComplexRateEigen {
    std::complex<double>* cc_ijk;       // This is a pre-computed tensor to avoid computing V D V^-1 every time. Cijk = Vik (V^-1)kj
    std::complex<double>* ceigenvalue;  // The complex eigenvalues obtained by EigenSystem
    std::complex<double>* oldCC_ijk;    // Old values to reset the current ones if there was a rejection
    std::complex<double>* oldCeigenvalue;
    int numStates;

    // Constructor
    ComplexRateEigen(int nS) : numStates(nS) {
        cc_ijk = new std::complex<double>[numStates * numStates * numStates];
        ceigenvalue = new std::complex<double>[numStates];
        oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
        oldCeigenvalue = new std::complex<double>[numStates];

        std::fill(cc_ijk, cc_ijk + numStates * numStates * numStates, 0.0);
        std::fill(oldCC_ijk, oldCC_ijk + numStates * numStates * numStates, 0.0);
        std::fill(ceigenvalue, ceigenvalue + numStates, 0.0);
        std::fill(oldCeigenvalue, oldCeigenvalue + numStates, 0.0);  
    }
    
    // Copy constructor
    ComplexRateEigen(const ComplexRateEigen& other) {
        numStates = other.numStates;
        cc_ijk = new std::complex<double>[numStates * numStates * numStates];
        ceigenvalue = new std::complex<double>[numStates];
        oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        
        std::copy(other.cc_ijk, other.cc_ijk + numStates * numStates * numStates, cc_ijk);
        std::copy(other.oldCC_ijk, other.oldCC_ijk + numStates * numStates * numStates, oldCC_ijk);
        std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
        std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);
    }

    // Assignment operator
    ComplexRateEigen& operator=(const ComplexRateEigen& other) {
        if (this != &other) {
            delete [] cc_ijk;
            delete [] ceigenvalue;
            delete [] oldCC_ijk;
            delete [] oldCeigenvalue;

            numStates = other.numStates;
            cc_ijk = new std::complex<double>[numStates * numStates * numStates];
            ceigenvalue = new std::complex<double>[numStates];
            oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
            oldCeigenvalue = new std::complex<double>[numStates];

            std::copy(other.cc_ijk, other.cc_ijk + numStates * numStates * numStates, cc_ijk);
            std::copy(other.oldCC_ijk, other.oldCC_ijk + numStates * numStates * numStates, oldCC_ijk);
            std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
            std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);
        }
        return *this;
    }

    ~ComplexRateEigen() {
        delete [] cc_ijk;
        delete [] ceigenvalue;
        delete [] oldCC_ijk;
        delete [] oldCeigenvalue;
    }
};

#endif