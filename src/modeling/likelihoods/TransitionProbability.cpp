#include <complex>
#include "TransitionProbability.hpp"
#include "core/Math.hpp"

/* This constructor builds a time-reversible transition matrix to describe
   evolutionary substitutions along a phylogenetic tree for discrete characters.
   Such a model is often used to describe substitutions for DNA or amino acid data
   but it works equally well for other characters with discrete states. */
TransitionProbability::TransitionProbability(const int nn, const int nC, const Matrix<double> &qMat)
    : ceigValExp(0), eigens(0), eigValExp(0), isComplex(false), isOldComplex(false),
	  numStates(qMat.dim1()), numNodes(nn), numCats(nC) {

	// Initialize
	numStates = qMat.dim1();
	initializeProbabilityBuffer();

	// Allocate space for variables that hold Q
	Q    = Matrix<double>(numStates, numStates, 0.0);
	pi = std::vector<double>(numStates, 0.0);

	ceigValExp = new std::complex<double>[numStates];
	eigens = new EigenSystem(Q);
	eigValExp = new double[numStates];
	
	updateQ(qMat, 0);
	accept();
}

/* Destructor. Deallocates memory used for Q matrix
   and eigensystem. */
TransitionProbability::~TransitionProbability(void) {

	delete [] ceigValExp;
	delete eigens;
	delete [] eigValExp;
	
	for(int i = 0; i < numNodes*numCats; i++) {
        delete probs[0][i];
        delete probs[1][i];
    }

    delete [] probs[0];
}

/* This function precalculates the product of the eigenvectors and their
   inverse for faster calculation of transition probabilities. The output
   is a vector of precalculated values (c_ijk). The input is the eigenvector
   matrix and the inverse of the eigenvector matrix. This function also
   fetches the eigenvalues from the eigensystem and stores them in an array
   of doubles in this class. */
void TransitionProbability::calcCijk(int mIndex) {

	// keep a copy of the eigenvalues
	double* p = &(eigens->getRealEigenvalues()[0]);
	double* q = rateEigen[mIndex].eigenvalue;
	memcpy(q, p, numStates*sizeof(double));

	// calculate c_ijk
	Matrix<double> ev  = eigens->getEigenvectors();
	Matrix<double> iev = eigens->getInverseEigenvectors();
	double* pc = rateEigen[mIndex].c_ijk;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			for (int k=0; k<numStates; k++)
			 	*(pc++) = ev(i, k) * iev(k, j);
}

/* This function precalculates the product of the eigenvectors and their
   inverse for faster calculation of transition probabilities when we have
   at least one complex eigenvalue. The output is a vector of precalculated
   complex values (cc_ijk). The input is the complex eigenvector matrix
   and the inverse of the complex eigenvector matrix. This function also
   fetches the real and imaginary eigenvalues from the eigensystem and stores
   them in two arrays of double values (eigenvalues and ieigenvalues) in this
   class. */
void TransitionProbability::calcComplexCijk(int mIndex) {

	// keep a copy of the complex eigenvalues
	double* p = &(eigens->getRealEigenvalues()[0]);
	double* q = &(eigens->getImagEigenvalues()[0]);
	for (int i=0; i<numStates; i++)
		complexRateEigen[mIndex].ceigenvalue[i] = std::complex<double>(*p++, *q++);

	// calculate cc_ijk
	Matrix<std::complex<double> > cev = eigens->getComplexEigenvectors();
	Matrix<std::complex<double> > ciev = eigens->getComplexInverseEigenvectors();
	std::complex<double>* pc = complexRateEigen[mIndex].cc_ijk;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			for (int k=0; k<numStates; k++)
			 	*(pc++) = cev(i, k) * ciev(k, j);
}

/* This function calculates the stationary frequencies of the rate matrix. The
   rate matrix, Q, is the infinitesimal generator of the Markov chain. It is an
   n X n matrix whose off-diagonal elements are q_ij >= 0 and whose diagonal elements
   are specified such that each row sums to zero. The rate matrix is finite (has
   a fixed number of states) and we assume that the input matrix is irreducible, as
   is the usual case for substitution models. Because Q is irreducible and finite,
   it has a stationary distribution, pi, which is a row vector of n probabilities.
   The stationary probabilities can be calculated by solving the homogeneous system
   of equations, pi*Q = 0, where 0 is a vector of zeros.
  
   We do the following to calculate the stationary frequencies.
  
   1. We perform an LU decomposition of the transpose of the matrix Q.
  
      Q' = LU
  
   2. Now we set Ux = z (x will eventually hold the stationary probabilities).
      Because L is nonsingular, we have z = 0. We proceed to back substitute on
      Ux = z = 0. When u_nn = 0, we can put in any solution for x. Here, we put
      in x_n = 1. We then solve the other values of x through back substitution.
  
   3. The solution obtained in 2 is not a probability vector. We normalize the
      vector such that the sum of the elements is 1.
  
   Note that the only time we need to use this function is when we don't
   know the stationary frequencies of the rate matrix beforehand. For most
   substitution models used in molecular evolution, the stationary frequencies
   are built into the rate matrix itself. These models are time-reversible.
   This function is useful for the non-reversible models.
  
   Stewart, W. J. 1999. Numerical methods for computing stationary distributions of
      finite irreducible Markov chains. In "Advances in Computational
      Probability", W. Grassmann, ed. Kluwer Academic Publishers. */
void TransitionProbability::calcStationaryFreq(void) {

	// transpose the rate matrix (qMatrix) and put into QT
	Matrix<double> QT(numStates, numStates);
	Math::transposeMatrix(Q, QT);

	// compute the LU decomposition of the transposed rate matrix
	Matrix<double> L(numStates, numStates);
	Matrix<double> U(numStates, numStates);
	Math::computeLandU(QT, L, U);
	
	// back substitute into z = 0 to find un-normalized stationary frequencies
	// start with x_n = 1.0
	pi[numStates-1] = 1.0;
	for (int i=numStates-2; i>=0; i--)
		{
		double dotProduct = 0.0;
		for (int j=i+1; j<numStates; j++)
			dotProduct += U(i, j) * pi[j];
		pi[i] = (0.0 - dotProduct) / U(i, i);
		}
		
	// normalize the solution vector
	double sum = 0.0;
	for (int i=0; i<numStates; i++)
		sum += pi[i];
	for (int i=0; i<numStates; i++)
		pi[i] /= sum;
}

void TransitionProbability::initializeProbabilityBuffer(){
	probs[0] = new Matrix<double>*[2*numNodes*numCats];
    probs[1] = probs[0] + numNodes;

    for(int i = 0; i < numNodes*numCats; i++){
        probs[0][i] = new Matrix<double>(numStates, numStates, 0.0);
        probs[1][i] = new Matrix<double>(numStates, numStates, 0.0);
    }
}

/* Rescale the rate matrix so that the mean rate at stationarity
   is 1.0. Requires that stationary frequencies have been
   calculated first. */
void TransitionProbability::rescaleQ (bool scalePi) {

	// Calculate the scaler, a factor by which all elements of Q
	// are multiplied such that the mean rate of substitution is 1
	double scaler = 0.0;
	for (int i=0; i<numStates; i++)
		if(scalePi)
			scaler += pi[i] * Q(i, i);
		else
			scaler += Q(i, i); // If the stationary is already included
			
	// Rescale rate matrix
	scaler = -1.0 / scaler;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			Q(i, j) *= scaler;
}

void TransitionProbability::accept(void) {
	isOldComplex = isComplex;
	for(int i = 0; i < isComplex.size(); i++){
		if(!isComplex[i]){
			memcpy(rateEigen[i].oldC_ijk, rateEigen[i].c_ijk, numStates*numStates*numStates*sizeof(double));
			memcpy(rateEigen[i].oldEigenvalue, rateEigen[i].eigenvalue, numStates*sizeof(double));
		}
		else {
			memcpy(complexRateEigen[i].oldCC_ijk, complexRateEigen[i].cc_ijk, numStates*numStates*numStates*sizeof(std::complex<double>));
			memcpy(complexRateEigen[i].oldCeigenvalue, complexRateEigen[i].ceigenvalue, numStates*sizeof(std::complex<double>));
		}
	}
}

void TransitionProbability::reject(void) {	
	isComplex = isOldComplex;
	for(int i = 0; i < isOldComplex.size(); i++){
		if(!isComplex[i]){
			memcpy(rateEigen[i].c_ijk, rateEigen[i].oldC_ijk, numStates*numStates*numStates*sizeof(double));
			memcpy(rateEigen[i].eigenvalue, rateEigen[i].oldEigenvalue, numStates*sizeof(double));
		}
		else {
			memcpy(complexRateEigen[i].cc_ijk, complexRateEigen[i].oldCC_ijk, numStates*numStates*numStates*sizeof(std::complex<double>));
			memcpy(complexRateEigen[i].ceigenvalue, complexRateEigen[i].oldCeigenvalue, numStates*sizeof(std::complex<double>));
		}
	}
}


/* This function returns transition probabilities in the Matrix P */
void TransitionProbability::setProbs(const int state, const int rate, const int node, const double v) {
	Matrix<double> P = *(probs[state][rate*numNodes + node]);
	if (!isComplex[rate])
		tiProbsEigens(v, P, rate);
	else
		tiProbsComplexEigens(v, P, rate);
}

/* This function calculates transition probabilities using
   complex eigenvalues and eigenvectors. */
void TransitionProbability::tiProbsComplexEigens(const double v, Matrix<double>& P, const int mIndex) {
	
	for (int s=0; s<numStates; s++)
		ceigValExp[s] = exp(complexRateEigen[mIndex].ceigenvalue[s] * v);

	const std::complex<double>* ptr = complexRateEigen[mIndex].cc_ijk;
	for (int i=0; i<numStates; i++)
		{
		for (int j=0; j<numStates; j++) 
			{
			std::complex<double> sum = std::complex<double>(0.0, 0.0);
			for(int s=0; s<numStates; s++)
				sum += (*ptr++) * ceigValExp[s];
			P(i, j) = (sum.real() < 0.0) ? 0.0 : sum.real();
			}
		}
}

/* This function calculates transition probabilities using
   eigenvalues and eigenvectors. */
void TransitionProbability::tiProbsEigens(const double v, Matrix<double> &P, const int mIndex) {
	
	for (int s=0; s<numStates; s++)
		eigValExp[s] = exp(rateEigen[mIndex].eigenvalue[s] * v);

	double *ptr = rateEigen[mIndex].c_ijk;
	for (int i=0; i<numStates; i++) 
		{
		for (int j=0; j<numStates; j++) 
			{
			double sum = 0.0;
			for(int s=0; s<numStates; s++)
				sum += (*ptr++) * eigValExp[s];
			P(i, j) = (sum < 0.0) ? 0.0 : sum;
			}
		}
}

int TransitionProbability::updateQ(const Matrix<double>& qTemp, const int index) {

	// Initialize the Q matrix
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			Q(i, j) = qTemp(i, j);
	
	rescaleQ(false);

	eigens->update(Q);

	// Do we need to make a new entry?
	if(index >= isComplex.size()) {
		isComplex.push_back(eigens->getIsComplex());
		rateEigen.push_back(RateEigen(numStates));
		complexRateEigen.push_back(ComplexRateEigen(numStates));
	}

	// Precalculate the product of the eigenvectors and their inverse
	if (isComplex[index] == false) {
			calcCijk(index);
	}
	else {
			calcComplexCijk(index);
	}

	return (0);
}

// Be sure you want to delete!!
void TransitionProbability::deleteQ(const int index) {
	isComplex.erase(isComplex.begin() + index);
	rateEigen.erase(rateEigen.begin() + index);
	complexRateEigen.erase(complexRateEigen.begin() + index);
}