#ifndef MATH_HPP
#define MATH_HPP

#include "Matrix.hpp"
#include <cmath>
#include <vector>

/**
 * @brief This namespace contains math utility functions. Call these functions
 * by using Math::<function>. An alternative is to declare 'using
 * namespace Math;', after which functions can be accessed without the
 * Math:: prefix. Some functions return 0 on success and 1 on failure;
 * it is up to the calling function to handle the error. 
 * 
 * @note All of this code was provided by John Huelsenbeck from a version of MrBayes.
 * The only exception is the choleskyDecomposition function, which was added to implement
 * Bayesian optimization.
 */
namespace Math {

    void    backSubstitutionRow(Matrix<double>& u, std::vector<double>& b);                         // Back substitution of row
    void    choleskyDecomposition(Matrix<double>& spd, Matrix<double>& cf);                         // Cholesky Decomposition
    void    computeLandU(Matrix<double>& aMat, Matrix<double>& lMat, Matrix<double>& uMat);         // LU decomposition
    int     expMatrixPade(Matrix<double>& a, Matrix<double>& f, int q);                             // Exponentiate matrix using Pade approximation
    double  factorial(int x);                                                                       // Return x! (x factorial)
    int     findPadeQValue(const double tolerance);                                                 // Set p and q of the Pade method to achieve desired tolerance
    void    forwardSubstitutionRow(Matrix<double>& L, std::vector<double>& b);                      // Forward substitution of row
    void    gaussianElimination (Matrix<double>& a, Matrix<double>& bMat, Matrix<double>& xMat);    // Gaussian elimination
    double  hypotenuse(double a, double b);                                                         // Return hypotenuse of triangle with legs a and b
    double  lnFactorial(int x);                                                                     // Return ln(x!)
    double  lnGamma(double alp);                                                                    // Return lnGamma(alp)
    int     transposeMatrix(const Matrix<double>& a, Matrix<double>& t);                            // Transpose matrix
    unsigned long stirlingFirst(int n, int k);
    unsigned long lnStirlingFirst(int n, int k);
}

#endif

