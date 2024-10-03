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
class EigenSystem {

	public:
                                EigenSystem(const Matrix<double>& m);                                                                                  //!< construct the eigenvalue decomposition
                               ~EigenSystem(void);                                                                                                       //!< destructor
        double                  getDeterminant(void);                                                                                                    //!< return determinant
        Matrix<double>&       getEigenvectors(void) { return eigenvectors; }                                                                           //!< return the eigenvector matrix
        Matrix<double>&       getInverseEigenvectors(void) { return inverseEigenvectors; }                                                             //!< return the inverse eigenvector matrix
        std::vector<double>&    getRealEigenvalues(void) { return realEigenvalues; }                                                                     //!< return the real parts of the eigenvalues
        std::vector<double>&    getImagEigenvalues(void) { return imaginaryEigenvalues; }                                                                //!< return the imaginary parts of the eigenvalues
        Matrix<complexNum>&      getComplexEigenvectors(void) { return complexEigenvectors; }                                                             //!< return the eigenvector matrix
        Matrix<complexNum>&      getComplexInverseEigenvectors() { return complexInverseEigenvectors; }                                                   //!< return the inverse eigenvector matrix
        bool                    getIsComplex(void) { return isComplex; }                                                                                 //!< returns 'true' if there are complex eigenvalues
        int                     update(const Matrix<double>& m);                                                                                       //!< update the eigensystem for matrix m

    private:
        int                     n;                                                                                                                                         //!< row and column dimension (square matrix)
        Matrix<double>        eigenvectors;                                                                                                            //!< matrix for internal storage of eigenvectors
        Matrix<double>        inverseEigenvectors;                                                                                                     //!< matrix for internal storage of the inverse eigenvectors
        Matrix<complexNum>       complexEigenvectors;                                                                                                     //!< matrix for internal storage of complex eigenvectors
        Matrix<complexNum>       complexInverseEigenvectors;                                                                                              //!< matrix for internal storage of the inverse of the complex eigenvectors
        std::vector<double>     realEigenvalues;                                                                                                         //!< vector for internal storage of the eigenvalues (real part)
        std::vector<double>     imaginaryEigenvalues;                                                                                                    //!< vector for internal storage of the eigenvalues (imaginary part)
        bool                    isComplex;                                                                                                               //!< flag whether there are complex eigenvalues
        void                    allocateComplexEigenvectors(void);                                                                                       //!< allocate space for complex eigenvectors
        void                    balance(Matrix<double>& A, std::vector<double>& scale, int* low, int* high);                                           //!< balances a matrix
        void                    balback(int low, int high, std::vector<double>& scale, Matrix<double>& eivec);                                         //!< reverses the balancing
        bool                    checkForComplexEigenvalues(void);                                                                                        //!< returns 'true' if there are complex eigenvalues
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

