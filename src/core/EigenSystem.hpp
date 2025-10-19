#ifndef EIGEN_SYSTEM_HPP
#define EIGEN_SYSTEM_HPP

#include "Matrix.hpp"
#include <complex>
#include <vector>

typedef std::complex<double> complexNum;

struct RateEigen;          // Defined in RateEigens.hpp
struct ComplexRateEigen;   // Defined in RateEigens.hpp

/**
 * @brief This class constructs the eigenvalue decomposition of
 * a square real (i.e., non-complex) matrix. The inverse of the
 *  eigenvectors are also calculated and stored.
 *
 *     Peters, G., and J.H. Wilkinson. 1970. Eigenvectors of real
 *        and complex matrices by LR and QR triangularisations.
 *        Numer. Math. 16:184-204.
 *     Martin, R.S., and J.H. Wilkinson. 1968. Similarity reduction
 *        of a general matrix to Hessenberg form. Numer. Math.
 *        12:349-368.
 *     Parlett, B.N., and C. Reinsch. 1969. Balancing a matrix for
 *        calculation of eigenvalues and eigenvectors. Numer.
 *        Math. 13:292-304. 
 * @note All of this code was provided by John Huelsenbeck from a version of MrBayes.
 * The only modifications that have been made are to allow for it to be compatible with
 * multi-threading (by no longer storing the eigen values) and use RateEigen/ComplexRateEigen.
*/
class EigenSystem {

	public:
                                EigenSystem(void) {};                                                                                                    // Construct the eigenvalue decomposition
                               ~EigenSystem(void);                                                                                                       // Destructor
        double                  getDeterminant(std::vector<double> realEigenValues);                                                                     // Return determinant
        bool                    update(const Matrix<double>& m, RateEigen& eigens, ComplexRateEigen& complexEigens);                                     // Update the eigensystem for matrix m

    private:                                                                                                                                             // Row and column dimension (square matrix)
        void                    allocateComplexEigenvectors(void);                                                                                       // Allocate space for complex eigenvectors
        void                    balance(Matrix<double>& A, std::vector<double>& scale, int* low, int* high);                                             // Balances a matrix
        void                    balback(int low, int high, std::vector<double>& scale, Matrix<double>& eivec);                                           // Reverses the balancing
        void                    complexLUBackSubstitution(Matrix<complexNum>& a, int *indx, std::vector<complexNum>& b);                                 // Back-substitutes a complex LU-decomposed matrix
        int                     complexLUDecompose(Matrix<complexNum>& a, double *vv, int* indx, double* pd);                                            // Calculates the LU-decomposition of a complex matrix
        void                    elmhes(int low, int high, Matrix<double>& a, std::vector<int>& perm);                                                    // Reduces matrix to upper Hessenberg form
        void                    elmtrans(int low, int high, Matrix<double> &a, std::vector<int> &perm, Matrix<double>& h);                               // Copies the Hessenberg matrix
        int                     hqr2(int low, int high, Matrix<double>& h, std::vector<double>& wr, std::vector<double>& wi, Matrix<double>& eivec);     // Computes eigenvalues and eigenvectors
        void                    initializeComplexEigenvectors(void);                                                                                     // Sets up the complex eigenvector matrix
        int                     invertMatrix(Matrix<double>& a, Matrix<double>& aInv);                                                                   // Inverts a matrix
        int                     invertComplexMatrix(Matrix<complexNum>& a, Matrix<complexNum>& aInv);                                                    // Inverts a complex matrix
        void                    luBackSubstitution (Matrix<double>& a, int* indx, double* b);                                                            // Back-substitutes an LU-decomposed matrix
        int                     luDecompose(Matrix<double>& a, double* vv, int* indx, double* pd);                                                       // Calculates the LU-decomposition of a matrix
};

#endif

