#ifndef RATE_EIGEN_HPP
#define RATE_EIGEN_HPP
#include <complex>
#include <vector>

struct RateEigen {
    double* c_ijk;
    double* eigenvalue;
    double* oldC_ijk;
    double* oldEigenvalue;
    int numStates;

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

struct ComplexRateEigen {
    std::complex<double>* cc_ijk;
    std::complex<double>* ceigenvalue;
    std::complex<double>* oldCC_ijk;
    std::complex<double>* oldCeigenvalue;
    int numStates;

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