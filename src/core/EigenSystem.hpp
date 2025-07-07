#ifndef EIGEN_SYSTEM_HPP
#define EIGEN_SYSTEM_HPP

#include "Matrix.hpp"
#include <complex>
#include <vector>

typedef std::complex<double> complexNum;

/* This class constructs and stores the eigenvalue decomposition of
   a square real (i.e., non-complex) matrix. The inverse of the
   eigenvectors are also calculated and stored.

      Peters, G., and J.H. Wilkinson. 1970. Eigenvectors of real
         and complex matrices by LR and QR triangularisations.
         Numer. Math. 16:184-204.
      Martin, R.S., and J.H. Wilkinson. 1968. Similarity reduction
         of a general matrix to Hessenberg form. Numer. Math.
         12:349-368.
      Parlett, B.N., and C. Reinsch. 1969. Balancing a matrix for
         calculation of eigenvalues and eigenvectors. Numer.
         Math. 13:292-304. */

struct RateEigen;
struct ComplexRateEigen;

class EigenSystem {

	public:
                                EigenSystem(void) {};                                                                                  //!< construct the eigenvalue decomposition
                               ~EigenSystem(void);                                                                                                       //!< destructor
        double                  getDeterminant(std::vector<double> realEigenValues);                                                                                                    //!< return determinant
        bool                    update(const Matrix<double>& m, RateEigen& eigens, ComplexRateEigen& complexEigens);                                                                                       //!< update the eigensystem for matrix m

    private:                                                                                                                                    //!< row and column dimension (square matrix)
        void                    allocateComplexEigenvectors(void);                                                                                       //!< allocate space for complex eigenvectors
        void                    balance(Matrix<double>& A, std::vector<double>& scale, int* low, int* high);                                           //!< balances a matrix
        void                    balback(int low, int high, std::vector<double>& scale, Matrix<double>& eivec);                                         //!< reverses the balancing
        void                    complexLUBackSubstitution(Matrix<complexNum>& a, int *indx, std::vector<complexNum>& b);                                     //!< back-substitutes a complex LU-decomposed matrix
        int                     complexLUDecompose(Matrix<complexNum>& a, double *vv, int* indx, double* pd);                                             //!< calculates the LU-decomposition of a complex matrix
        void                    elmhes(int low, int high, Matrix<double>& a, std::vector<int>& perm);                                                  //!< reduces matrix to upper Hessenberg form
        void                    elmtrans(int low, int high, Matrix<double> &a, std::vector<int> &perm, Matrix<double>& h);                           //!< copies the Hessenberg matrix
        int                     hqr2(int low, int high, Matrix<double>& h, std::vector<double>& wr, std::vector<double>& wi, Matrix<double>& eivec); //!< computes eigenvalues and eigenvectors
        void                    initializeComplexEigenvectors(void);                                                                                     //!< sets up the complex eigenvector matrix
        int                     invertMatrix(Matrix<double>& a, Matrix<double>& aInv);                                                               //!< inverts a matrix
        int                     invertComplexMatrix(Matrix<complexNum>& a, Matrix<complexNum>& aInv);                                                      //!< inverts a complex matrix
        void                    luBackSubstitution (Matrix<double>& a, int* indx, double* b);                                                          //!< back-substitutes an LU-decomposed matrix
        int                     luDecompose(Matrix<double>& a, double* vv, int* indx, double* pd);                                                     //!< calculates the LU-decomposition of a matrix
};

#endif

