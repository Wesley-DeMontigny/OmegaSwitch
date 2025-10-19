#ifndef MATRIX_HPP 
#define MATRIX_HPP 
 
#include <iostream> 
#include <iomanip> 
#include <cstdlib> 
#include <cstring>
 
// We need to forward declare some stuff
template <class T>
class Matrix;

template <class T> std::ostream& 	operator<<(std::ostream &s, const Matrix<T> &A);
template <class T> std::istream& 	operator>>(std::istream &s, Matrix<T> &A);
template <class T> Matrix<T> 		operator+(const Matrix<T> &A, const Matrix<T> &B);
template <class T> Matrix<T> 		operator-(const Matrix<T> &A, const Matrix<T> &B);
template <class T> Matrix<T> 		operator*(const Matrix<T> &A, const Matrix<T> &B);
template <class T> Matrix<T> 		&operator+=(Matrix<T> &A, const Matrix<T> &B);
template <class T> Matrix<T> 		&operator-=(Matrix<T> &A, const Matrix<T> &B);
template <class T> Matrix<T> 		&operator*=(Matrix<T> &A, const Matrix<T> &B);
template <class T> Matrix<T> 		operator+(const T &a, const Matrix<T> &B);
template <class T> Matrix<T> 		operator-(const T &a, const Matrix<T> &B);
template <class T> Matrix<T> 		operator*(const T &a, const Matrix<T> &B);
template <class T> Matrix<T> 		operator/(const T &a, const Matrix<T> &B);
template <class T> Matrix<T> 		operator+(const Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		operator-(const Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		operator*(const Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		operator/(const Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		&operator+=(Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		&operator-=(Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		&operator*=(Matrix<T> &A, const T &b);
template <class T> Matrix<T> 		&operator/=(Matrix<T> &A, const T &b);


/**
 * @brief This class implements all of the constructors, destructors, and operators required
 * for us to do basic Matrix/vector algebra
 * 
 * @tparam T The type the matrix is made up of
 * @note All of this code was provided by John Huelsenbeck from a version of MrBayes. The only
 * modifications that have been done to it are that the matrix is now accessed using the () operators,
 * which can take both a row index and a column index, as opposed to the previous [] operator.
 */
template <class T> 
class Matrix { 
 
public: 
		Matrix(void);                                   //!< null constructor (0 X 0 matrix) 
		Matrix(int m, int n);                           //!< creates a m X n matrix without initialization 
		Matrix(int m, int n, const T &a);               //!< creates a m X n matrix initializing all elements to a constant 
		inline Matrix(const Matrix &A);                      //!< creates a m X n matrix with elements shared by another matrix, A 
		~Matrix(void);                                  //!< destructor 
		operator T*() { return v; }                //!< type cast operator
		operator const T*() { return v; }          //!< type cast operator for const
		inline Matrix   &operator=(const T &a);                           //!< assignment operator (all elements have the value a) 
		Matrix   		&operator=(const Matrix &A) { return ref(A); }  //!< assignment operator (shallow copy, elements share data) 
		bool   			operator==(const Matrix &A) const;              //!< equality operator 
		T&   			operator()(size_t r, size_t c) { return this->v[r*n + c]; }
		const T&   		operator()(size_t r, size_t c) const { return this->v[r*n + c]; }
		friend std::ostream& 	operator<<<>(std::ostream &s, const Matrix<T> &A);  //!< operator << 
		friend std::istream& 	operator>><>(std::istream &s, Matrix<T> &A);        //!< operator >> 
		friend Matrix<T> 		operator+<>(const Matrix<T> &A, const Matrix<T> &B);     //!< operator + 
		friend Matrix<T> 		operator-<>(const Matrix<T> &A, const Matrix<T> &B);     //!< operator - 
		friend Matrix<T> 		operator*<>(const Matrix<T> &A, const Matrix<T> &B);     //!< operator * (matrix multiplication) 
		friend Matrix<T> 		&operator+=<>(Matrix<T> &A, const Matrix<T> &B);   //!< operator += 
		friend Matrix<T> 		&operator-=<>(Matrix<T> &A, const Matrix<T> &B);   //!< operator -= 
		friend Matrix<T> 		&operator*=<>(Matrix<T> &A, const Matrix<T> &B);   //!< operator *= (matrix multiplication)
		friend Matrix<T> 		operator+<>(const T &a, const Matrix<T> &B);               //!< operator + for scalar + matrix 
		friend Matrix<T> 		operator-<>(const T &a, const Matrix<T> &B);               //!< operator - for scalar - matrix 
		friend Matrix<T> 		operator*<>(const T &a, const Matrix<T> &B);               //!< operator * for scalar * matrix 
		friend Matrix<T> 		operator/<>(const T &a, const Matrix<T> &B);               //!< operator / for scalar / matrix 
		friend Matrix<T> 		operator+<>(const Matrix<T> &A, const T &b);               //!< operator + for matrix + scalar 
		friend Matrix<T> 		operator-<>(const Matrix<T> &A, const T &b);               //!< operator - for matrix - scalar 
		friend Matrix<T> 		operator*<>(const Matrix<T> &A, const T &b);               //!< operator * for matrix * scalar 
		friend Matrix<T> 		operator/<>(const Matrix<T> &A, const T &b);               //!< operator / for matrix / scalar 
		friend Matrix<T> 		&operator+=<>(Matrix<T> &A, const T &b);             //!< operator += for scalar 
		friend Matrix<T> 		&operator-=<>(Matrix<T> &A, const T &b);             //!< operator -= for scalar 
		friend Matrix<T> 		&operator*=<>(Matrix<T> &A, const T &b);             //!< operator *= for scalar 
		friend Matrix<T> 		&operator/=<>(Matrix<T> &A, const T &b);             //!< operator /= for scalar 
		inline Matrix   &ref(const Matrix &A);                          //!< creates a reference to another matrix (shallow copy) 
		Matrix   		copy(void) const;                                 //!< creates a copy of another matrix (deep copy, with separate data elements) 
		Matrix   		&inject(const Matrix &A);                       //!< copy the values of elements from one matrix to another 
		int   			dim1(void) const { return m; }                    //!< number of rows 
		int   			dim2(void) const { return n; }                    //!< number of columns 
		T*   			expose(void) { return v; }
		int   			getRefCount(void) const { return *refCount; }     //!< get the number of matrices that share the same data 

	private: 
		T*				v; 
		int   			m; 
		int   			n; 
		int*			refCount; 
		bool   			cArray; 

		void   			destroy(void); 
 
}; 

template <class T> bool        operator!=(const Matrix<T> &A, const Matrix<T> &B);   //!< inequality (operator !=) 

// Definitions of inlined member functions

/*!
 * Copy constructor, which creates a shallow copy of the
 * MbMatrix argument. Matrix data are not copied but shared.
 * Thus, in MbMatrix B(A), subsequent changes to A will be
 * reflected by changes in B. For an independent copy, use
 * MbMatrix B(A.copy()), or B = A.copy(), instead. Note
 * the use of garbage collection in this class, through the
 * reference counter refCount.
 *
 * \brief Shallow copy constructor
 * \param A Matrix to copy
 */
template <class T>
inline Matrix<T>::Matrix(const Matrix &A)
    : v(A.v), m(A.m), n(A.n), refCount(A.refCount), cArray(A.cArray) {

	(*refCount)++;
}

/*!
 * Assign all elements of the matrix the value of
 * the constant scalar a.
 *
 * \brief Assign scalar to all elements
 * \param a Scalar used in assignment
 * \return Assigned matrix
 */
template <class T>
inline Matrix<T> &Matrix<T>::operator=(const T &a) {

	for (int i =0, len = m*n; i<len; i++) {
		v[i] = a;
	}
	return *this;

}

/*!
 * Create a reference (shallow assignment) to another existing matrix.
 * In B.ref(A), B and A share the same data and subsequent changes
 * to the matrix elements of one will be reflected in the other. Note that
 * the reference counter is always allocated, even for null matrices,
 * so we need not test whether refCount is NULL.
 *
 * This is what operator= calls, and B=A and B.ref(A) are equivalent
 * operations.
 *
 * \brief Make this reference to A
 * \param A Matrix to take reference of
 * \return Matrix to which the reference was assigned
 */
template <class T>
inline Matrix<T> &Matrix<T>::ref(const Matrix &A) {

	if (this != &A)
		{
		(*refCount)--;
		if (*refCount < 1)
			destroy();
		m = A.m;
		n = A.n;
		v = A.v;
		refCount = A.refCount;
		cArray = A.cArray;
        (*refCount)++;
		}
	return *this;
	
}


// Defintions of member functions that are not inlined

/*!
 * Null constructor. Creates a (0 X O) ('null') matrix.
 * Note that the reference count will be set to 1 for this
 * null matrix. This is to simplify the rest of the
 * code at the cost of allocating and deleting an int
 * everytime a null matrix is needed.
 *
 * \brief Null constructor
 */
template <class T>
Matrix<T>::Matrix(void)
    : v(0), m(0), n(0), refCount(0), cArray(false) {

	refCount = new int;
	*refCount = 1;
	
}

/*!
 * Create a new (m X n) matrix, without initializing matrix 
 * elements. If m or n are not positive, a null matrix is created. Note
 * that the reference count will be set to 1 regardless of whether
 * we create a null matrix. This is to simplify the rest of the
 * code at the cost of allocating and deleting an int every time
 * a null matrix is needed.
 *
 * This version avoids the O(m*n) initialization overhead.
 *
 * \brief Constructor of uninitialized matrix
 * \param m The first (row) dimension of the matrix
 * \param n The second (column) dimension of the matrix
 */
template <class T>
Matrix<T>::Matrix(int m, int n)
    : v(0), m(0), n(0), refCount(0), cArray(false) {

	if (m > 0 && n > 0) {
		// allocate and initialize pointers
		v = new T[m * n];
		this->m = m;
		this->n = n;
	}    
	refCount = new int;
	*refCount = 1;

}

/*!
 * Create a new (m X n) matrix, initializing matrix elements
 * to the constant value specified by the third argument. Most
 * often used to create a matrix of zeros, as in MbMatrix A(m, n, 0.0).
 * 
 * \brief Constructor of initialized matrix.
 * \param m The first (row) dimension of the matrix
 * \param n The second (column) dimension of the matrix
 * \param a Value for initialization.
 */
template <class T>
Matrix<T>::Matrix(int m, int n, const T &a)
    : v(0), m(0), n(0), refCount(0), cArray(false) {

	if (m > 0 && n > 0) {
		// allocate and initialize pointers
		v = new T[m*n];
		for (int i =0, len = m*n; i<len; i++) {
			v[i] = a;
		}
		this->m = m;
		this->n = n;
	}

	refCount = new int;
	*refCount = 1;
}

/*!
 * Destructor. Note that refCount is decreased and only if
 * refCount reaches 0 do we delete allocated memory. This is
 * garbage collection as implemented in JAVA and other languages.
 * Note that null matrices also have a reference count allocated,
 * so that we can always access the value of refCount.
 *
 * \brief Destructor with garbage collection
 */
template <class T>
Matrix<T>::~Matrix(void) {

	(*refCount)--;
	if (*refCount < 1)
		destroy();

}

/*!
 * Equality operator. The dimensions of the two matrices
 * are first compared. If they are not the same, then false
 * is returned. Second, all elements are compared. If
 * they are the same, true is returned, otherwise false
 * is returned. Note that this operator is not useful
 * for float and double matrices, but it is handy for int
 * and bool matrices, as well as for matrices of other types
 * that have a sensible operator!= defined.
 *
 * \brief Equality operator
 * \param A Matrix to compare (*this) to
 * \return True if (*this)==A, false otherwise.
 */
template <class T>
bool Matrix<T>::operator==(const Matrix &A) const {

	if (m != A.m || n != A.n)
		return false;
	T* p1 = v;
	T* p2 = A.v;
	T* end = p1+(m*n);
	for (; p1<end; p1++, p2++) {
		if (*p1 != *p2)
			return false;
	}
	return true;

}

/*!
 * Create a new version of an existing matrix.  Used in B = A.copy()
 * or in the construction of a new matrix that does not share
 * data with the copied matrix, e.g. in MbMatrix B(A.copy()).
 *
 * \brief Create independent copy
 * \return Copy of this.
 */
template <class T>
Matrix<T> Matrix<T>::copy(void) const {
	Matrix A(m, n);
	memcpy (A.v, v, m*n*sizeof(T));
	return A;
}

/*
 * Copy the elements from one matrix to another, in place.
 * That is, if you call B.inject(A), both A and B must conform
 * (i.e. have the same row and column dimensions).
 *
 * This differs from B = A.copy() in that references to B
 * before this assignment are also affected.  That is, if
 * we have 
 *
 * MbMatrix A(n);
 * MbMatrix C(n);
 * MbMatrix B(C);        (elements of B and C are shared) 
 *
 * then B.inject(A) affects both B and C, while B=A.copy() creates
 * a new matrix B which shares no data with C or A.
 *
 * A is the matrix from which elements will be copied.
 * The function returns an instance of the modifed matrix. That is, in 
 * B.inject(A), it returns B.  If A and B are not conformant, no 
 * modifications to B are made.
 *
 * \brief Inject elements of A into (*this)
 * \param A Matrix with elements to inject
 * \return Injected matrix
 */
template <class T>
Matrix<T> &Matrix<T>::inject(const Matrix &A) {

	if (A.m == m && A.n == n)
		memcpy(v, A.v, m*n*sizeof(T));
	return *this;
}

/*!
 * This is a garbage collector, which is
 * called only when there is no more element
 * referencing this matrix.
 *
 * \brief Garbage collection
 */
template <class T>
void Matrix<T>::destroy(void) {

	if (v != 0) {
		delete [] v;
	}
	if (refCount != 0)
		delete refCount;
	cArray = false;
}


// Definitions of related templated functions on matrices

/*!
 * Printing of a matrix to an ostream object.
 * We use the format:
 * [<m>,<n>]
 * ((v_11,v_12,v_13,...,v_1n),
 * (v_i1,v_i2,v_i3,...,v_in),
 * (v_m1,v_m2,v_m3,...,v_mn)) 
 *
 * \brief operator<<
 * \param A Matrix to output
 * \param s ostream to output to
 * \return ostream object (for additional printing)
 */
template <class T>
std::ostream& operator<<(std::ostream &s, const Matrix<T> &A) {

	int M = A.n;
	int N = A.m;
	s << "[" << M << "," << N << "]\n";
	s << "(";
	for (int i=0; i<M; i++) {
		s << "(";
		for (int j=0; j<N; j++) {
			s << A(i, j);
			if (j != N-1)
				s << ",";
		}
		if (i != M-1)
			s << "),\n";
		else
			s << ")";
	}
	s << ")";
	return s;
	
}

/*!
 * Reading of a matrix from an istream object.
 * We expect the format:
 * [<m>,<n>]
 * ((v_11,v_12,v_13,...,v_1n),
 * (v_i1,v_i2,v_i3,...,v_in),
 * (v_m1,v_m2,v_m3,...,v_mn)) 
 * On failure, a null matrix is returned.
 *
 * \brief operator>>
 * \param A Matrix to receive input
 * \param s istream to read from
 * \return istream object (for additional reading)
 */
template <class T>
std::istream& operator>>(std::istream &s, Matrix<T> &A) {

	A = Matrix<T>();	// make sure we return null matrix on failure
	int M, N;
	char c;
	s >> c;
	if (c != '[')
		return s;
	s >> M;
	s >> c;
	if (c != ',')
		return s;
	s >> N;
	Matrix<T> B(M,N);
	s >> c;
	if (c != ']')
		return s;
	s.ignore();  // ignore the newline
	s >> c;
	if (c != '(')
		return s;
	for (int i=0; i<M; i++) {
		s >> c;
		if (c != '(')
			return s;
		for (int j=0; j<N; j++) {
			s >> B(i, j);
			if (j < N-1) {
				s >> c;
				if (c != ',') return s;
			}
		}
		s >> c;
		if (c != ')')
			return s;
		if (i != M-1)
			s.ignore();  // ignore newline
	}
	s >> c;
	if (c != ')')
		return s;
	A = B;
	return s;

}

/*!
 * This function performs addition of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator+ (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A + b
 */
template <class T>
Matrix<T> operator+(const Matrix<T> &A, const T &b) {

	Matrix<T> B(A.copy());
	T* p1 = B.v;
	T* p2 = A.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = *p2 + b;
	return B;
}

/*!
 * This function performs subtraction of a scalar from
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator- (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A - b
 */
template <class T>
Matrix<T> operator-(const Matrix<T> &A, const T &b) {

	Matrix<T> B(A.copy());
	T* p1 = B.v;
	T* p2 = A.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = *p2 - b;
	return B;

}

/*!
 * This function performs multiplication of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator* (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A * b
 */
template <class T>
Matrix<T> operator*(const Matrix<T> &A, const T &b) {

	Matrix<T> B(A.copy());
	T* p1 = B.v;
	T* p2 = A.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = *p2 * b;
	return B;

}

/*!
 * This function performs division with a scalar of
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator/ (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A / b
 */
template <class T>
Matrix<T> operator/(const Matrix<T> &A, const T &b) {

	Matrix<T> B(A.copy());
	T* p1 = B.v;
	T* p2 = A.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = *p2 / b;
	return B;

}

/*!
 * This function performs addition of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator+ (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a + B
 */
template <class T>
Matrix<T> operator+(const T &a, const Matrix<T> &B) {

	Matrix<T> A(B.copy());
	T* p1 = A.v;
	T* p2 = B.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = a + *p2;
	return A;

}

/*!
 * This function subtracts each element of a
 * a matrix from a scalar and returns the
 * resulting matrix.
 *
 * \brief operator- (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a - B
 */
template <class T>
Matrix<T> operator-(const T &a, const Matrix<T> &B) {

	Matrix<T> A(B.copy());
	T* p1 = A.v;
	T* p2 = B.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = a - *p2;
	return A;

}

/*!
 * This function performs multiplication of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator* (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a * B
 */
template <class T>
Matrix<T> operator*(const T &a, const Matrix<T> &B) {

	Matrix<T> A(B.copy());
	T* p1 = A.v;
	T* p2 = B.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = a **p2;
	return A;

}

/*!
 * This function performs division of a scalar by
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator/ (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a / B
 */
template <class T>
Matrix<T> operator/(const T &a, const Matrix<T> &B) {

	Matrix<T> A(B.copy());
	T* p1 = A.v;
	T* p2 = B.v;
	T* end = p1 + (B.m*B.n);
	for (; p1<end; p1++, p2++)
		*p1 = a / *p2;
	return A;

}

/*!
 * This function performs addition of a scalar to
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator+= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A += b
 */
template <class T>
Matrix<T> &operator+=(Matrix<T> &A, const T &b) {

	T* p1 = A.v;
	T* end = p1 + (A.m*A.n);
	for (; p1<end; p1++)
		*p1 += b;
	return A;

}

/*!
 * This function performs subtraction of a scalar from
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator-= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A -= b
 */
template <class T>
Matrix<T> &operator-=(Matrix<T> &A, const T &b) {

	T* p1 = A.v;
	T* end = p1 + (A.m*A.n);
	for (; p1<end; p1++)
		*p1 -= b;
	return A;

}

/*!
 * This function performs multiplication of a scalar to
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator*= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A *= b
 */
template <class T>
Matrix<T> &operator*=(Matrix<T> &A, const T &b) {

	T* p1 = A.v;
	T* end = p1 + (A.m*A.n);
	for (; p1<end; p1++)
		*p1 *= b;
	return A;

}

/*!
 * This function performs division with a scalar of
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator/= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A /= b
 */
template <class T>
Matrix<T> &operator/=(Matrix<T> &A, const T &b) {

	T* p1 = A.v;
	T* end = p1 + (A.m*A.n);
	for (; p1<end; p1++)
		*p1 /= b;
	return A;
}

/*!
 * This function performs elementwise addition of two
 * matrices and returns the resulting matrix. If the
 * matrices are not conformant, a null matrix is returned.
 *
 * \brief operator+
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A + B, null matrix on failure
 */
template <class T>
Matrix<T> operator+(const Matrix<T> &A, const Matrix<T> &B) {

	int m = A.m;
	int n = A.n;
	if (B.m != m ||  B.n != n)
		return Matrix<T>();
	else {
		Matrix<T> C(m,n);
		T* p1 = A.v;
		T* p2 = B.v;
		T* p3 = C.v;
		T* end = p3 + (C.m*C.n);
		for (; p3<end; p1++, p2++, p3++)
			*p3 = *p1 + *p2;
		return C;
	}

}

/*!
 * This function performs elementwise subtraction of two
 * matrices and returns the resulting matrix. If the
 * matrices are not conformant, a null matrix is returned.
 *
 * \brief operator-
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A - B, null matrix on failure
 */
template <class T>
Matrix<T> operator-(const Matrix<T> &A, const Matrix<T> &B) {

	int m = A.m;
	int n = A.n;
	if (B.m != m ||  B.n != n)
		return Matrix<T>();
	else {
		Matrix<T> C(m,n);
		T* p1 = A.v;
		T* p2 = B.v;
		T* p3 = C.v;
		T* end = p3 + (C.m*C.n);
		for (; p3<end; p1++, p2++, p3++)
			*p3 = *p1 - *p2;
		return C;
	}

}

/*!
 * Compute C = A*B, where C[i][j] is the dot-product of 
 * row i of A and column j of B. Note that this operator
 * does not perform elementwise multiplication. If the 
 * matrices do not have the right dimensions for matrix
 * multiplication (that is, if the number of columns of A
 * is different from the number of rows of B), the function
 * returns a null matrix.
 *
 * \brief Matrix multiplication
 * \param A An (m X n) matrix
 * \param B An (n X k) matrix
 * \return A * B, an (m X k) matrix, or null matrix on failure
 */
template <class T>
Matrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B) {

	if ( A.n != B.m)
		return Matrix<T>();
	int M = A.m;
	int N = A.n;
	int K = B.n;
	Matrix<T> C(M,K);
	for (int i=0; i<M; i++) {
		for (int j=0; j<K; j++) {
			T sum = 0;
			for (int k=0; k<N; k++)
				sum += A(i, k) * B (k, j);
			C(i, j) = sum;
		}
	}
	return C;

}

/*!
 * This function performs elementwise addition on two
 * matrices and puts the result in the first matrix.
 * If the two matrices are nonconformant, the first
 * matrix is left intact.
 *
 * \brief operator+=
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A += B, A unmodified on failure
 */
template <class T>
Matrix<T>&  operator+=(Matrix<T> &A, const Matrix<T> &B) {

	int m = A.m;
	int n = A.n;
	if (B.m == m && B.n == n) {
		T* p1 = A.v;
		T* p2 = B.v;
		T* end = p1 + (A.m*A.n);
		for (; p1<end; p1++, p2++)
			*p1 += *p2;
	}
	return A;

}

/*!
 * This function performs elementwise subtraction on two
 * matrices and puts the result in the first matrix.
 * If the two matrices are nonconformant, the first
 * matrix is left intact.
 *
 * \brief operator-=
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A -= B, A unmodified on failure
 */
template <class T>
Matrix<T>&  operator-=(Matrix<T> &A, const Matrix<T> &B) {

	int m = A.m;
	int n = A.n;
	if (B.m == m && B.n == n) {
		T* p1 = A.v;
		T* p2 = B.v;
		T* end = p1 + (A.m*A.n);
		for (; p1<end; p1++, p2++)
			*p1 -= *p2;
	}
	return A;

}

/*!
 * Compute C = A*B, where C[i][j] is the dot-product of 
 * row i of A and column j of B. Then assign the result to
 * A. Note that this operator does not perform elementwise
 * multiplication. If the matrices are not both square of the
 * same dimension, then the operation is not possible to
 * perform and we return an unomidified A.
 *
 * \brief Matrix multiplication with assignment (operator *=)
 * \param A An (n X n) matrix
 * \param B An (n X n) matrix
 * \return A = A * B, an (n X n) matrix, or unmodified A on failure
 */
template <class T>
Matrix<T> &operator*=(Matrix<T> &A, const Matrix<T> &B) {

	if (A.m==A.n && B.m==B.n && A.m==B.m) {
		int N = A.m;
		Matrix<T> C(N,N);
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				T sum = 0;
				for (int k=0; k<N; k++)
					sum += A(i, k) * B (k, j);
				C(i, j) = sum;
			}
		}
		A = C;
	}
	return A;

}

/*!
 * Inequality operator. It calls operator== and returns
 * the reverse of the bool result. Note that this operator
 * is not useful for float and double matrices, but it is handy
 * for int and bool matrices, as well as for matrices of other types
 * that have a sensible operator!= defined.
 *
 * \brief Inequality operator
 * \param A Matrix 1
 * \param B Matrix 2
 * \return True if A != B, false otherwise
 */
template <class T>
bool operator!=(const Matrix<T> &A, const Matrix<T> &B) {

	if (A == B)
		return false;
	else
		return true;

}

#endif
