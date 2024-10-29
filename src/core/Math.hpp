#ifndef MATH_HPP
#define MATH_HPP

#include "Matrix.hpp"
#include <cmath>
#include <vector>

/*!
 * This namespace contains math utility functions. Call these functions
 * by using MbMath::<function>. An alternative is to declare 'using
 * namespace MbMath;', after which functions can be accessed without the
 * MbMath:: prefix. Some functions return 0 on success and 1 on failure;
 * it is up to the calling function to handle the error. */
namespace Math {

    void    backSubstitutionRow(Matrix<double>& u, std::vector<double>& b);                            //!< back substitution of row
    void    computeLandU(Matrix<double>& aMat, Matrix<double>& lMat, Matrix<double>& uMat);        //!< LU decomposition
    int     expMatrixPade(Matrix<double>& a, Matrix<double>& f, int q);                              //!< exponentiate matrix using Pade approximation
    double  factorial(int x);                                                                            //!< return x! (x factorial)
    int     findPadeQValue(const double tolerance);                                                      //!< set p and q of the Pade method to achieve desired tolerance
    void    forwardSubstitutionRow(Matrix<double>& L, std::vector<double>& b);                         //!< forward substitution of row
    void    gaussianElimination (Matrix<double>& a, Matrix<double>& bMat, Matrix<double>& xMat);   //!< gaussian elimination
    double  hypotenuse(double a, double b);                                                              //!< return hypotenuse of triangle with legs a and b
    double  lnFactorial(int x);                                                                          //!< return ln(x!)
    double  lnGamma(double alp);                                                                         //!< return lnGamma(alp)
    int     transposeMatrix(const Matrix<double>& a, Matrix<double>& t);                             //!< transpose matrix
    unsigned long stirlingFirst(int n, int k);
    unsigned long lnStirlingFirst(int n, int k);
}

#endif

