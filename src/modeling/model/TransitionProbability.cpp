#include <complex>
#include "TransitionProbability.hpp"
#include "core/Math.hpp"

TransitionProbability::TransitionProbability(const int nn, const int nC)
    : ceigValExp(0), eigens(0), eigValExp(0), isComplex(false), isOldComplex(false),
	  numStates(122), numNodes(nn), numCats(nC) {

	probs[0] = new Matrix<double>*[2*numNodes*numCats];
    probs[1] = probs[0] + numNodes;

    for(int i = 0; i < numNodes*numCats; i++){
        probs[0][i] = new Matrix<double>(numStates, numStates, 0.0);
        probs[1][i] = new Matrix<double>(numStates, numStates, 0.0);
    }

	Matrix<double> Q(122, 122, 0.0);

	ceigValExp = new std::complex<double>[numStates];
	eigens = new EigenSystem(Q);
	eigValExp = new double[numStates];
	
	updateQ(Q, 0);
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
	Matrix<std::complex<double>> cev = eigens->getComplexEigenvectors();
	Matrix<std::complex<double>> ciev = eigens->getComplexInverseEigenvectors();
	std::complex<double>* pc = complexRateEigen[mIndex].cc_ijk;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			for (int k=0; k<numStates; k++)
			 	*(pc++) = cev(i, k) * ciev(k, j);
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

void TransitionProbability::setProbs(const int state, const int rate, const int node, const double v) {
	Matrix<double> P = *(probs[state][rate*numNodes + node]);
	if (!isComplex[rate])
		tiProbsEigens(v, P, rate);
	else
		tiProbsComplexEigens(v, P, rate);
}

void TransitionProbability::pullProbs(const int state, const int rate, const int node, const double v) {
	Matrix<double> P = *(probs[state][rate*numNodes + node]);
	Matrix<double> P2 = *(probs[state ^ true][rate*numNodes + node]);

	for (int i=0; i<numStates; i++) {
		for (int j=0; j<numStates; j++) {
			P(i,j) = P2(i,j);
		}
	}
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

void TransitionProbability::updateQ(Matrix<double>& Q, const int index) {
	
	double scaler = 0.0;
	for (int i=0; i<numStates; i++)
		scaler += Q(i, i);
			
	// Rescale rate matrix
	scaler = -1.0 / scaler;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			Q(i, j) *= scaler;

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
}

// Be sure you want to delete!!
void TransitionProbability::deleteQ(const int index) {
	isComplex.erase(isComplex.begin() + index);
	rateEigen.erase(rateEigen.begin() + index);
	complexRateEigen.erase(complexRateEigen.begin() + index);
}