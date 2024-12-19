#include <complex>
#include "TransitionProbability.hpp"
#include "core/Math.hpp"
#include <cstring>

TransitionProbability::TransitionProbability(const int nn, const int nC)
    : numStates(122), numNodes(nn), numCats(nC) {

	probs[0] = new Matrix<double>*[2*numNodes*numCats];
    probs[1] = probs[0] + numNodes;

    for(int i = 0; i < numNodes*numCats; i++){
        probs[0][i] = new Matrix<double>(numStates, numStates, 0.0);
        probs[1][i] = new Matrix<double>(numStates, numStates, 0.0);
    }

	Matrix<double> Q(122, 122, 0.0);

	eigens = new EigenSystem(122);
	
	allocateQ(1);
	updateQ(Q, 0);
	accept();
}

/* Destructor. Deallocates memory used for Q matrix
   and eigensystem. */
TransitionProbability::~TransitionProbability(void) {

	delete eigens;
	
	for(int i = 0; i < numNodes*numCats; i++) {
        delete probs[0][i];
        delete probs[1][i];
    }

    delete [] probs[0];
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

	std::vector<std::complex<double>> ceigValExp;

	for (int s=0; s<numStates; s++)
		ceigValExp.push_back(exp(complexRateEigen[mIndex].ceigenvalue[s] * v));

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
	
	std::vector<double> eigValExp;

	for (int s=0; s<numStates; s++)
		eigValExp.push_back(exp(rateEigen[mIndex].eigenvalue[s] * v));

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

void TransitionProbability::allocateQ(int size){
	if(size > isComplex.size()) {
		for(int i = 0, num = size - isComplex.size(); i < num; i++){
			isComplex.push_back(false);
			rateEigen.push_back(RateEigen(numStates));
			complexRateEigen.push_back(ComplexRateEigen(numStates));
		}
	}
}

void TransitionProbability::updateQ(Matrix<double> Q, const int index) {
	
	double scaler = 0.0;
	for (int i=0; i<numStates; i++)
		scaler += Q(i, i);
			
	// Rescale rate matrix
	scaler = -1.0 / scaler;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			Q(i, j) *= scaler;

	isComplex[index] = eigens->update(Q, rateEigen[index], complexRateEigen[index]);
}

std::vector<Matrix<double>> TransitionProbability::generateProbs(Matrix<double> Q, std::vector<double> branches) {
	
	double scaler = 0.0;
	for (int i=0; i<numStates; i++)
		scaler += Q(i, i);
			
	// Rescale rate matrix
	scaler = -1.0 / scaler;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			Q(i, j) *= scaler;

	RateEigen localRateEigen(numStates);
	ComplexRateEigen localComplexRateEigen(numStates);

	bool isComplex = eigens->update(Q, localRateEigen, localComplexRateEigen);

	std::vector<Matrix<double>> returnMatrices;

	if(!isComplex) {
		for(int i = 0; i < branches.size(); i++){
			double v = branches[i];
			std::vector<double> eigValExp;

			Matrix<double> branchMatrix(numStates, numStates); 

			for (int s=0; s<numStates; s++)
				eigValExp.push_back(exp(localRateEigen.eigenvalue[s] * v));

			double *ptr = localRateEigen.c_ijk;
			for (int i=0; i<numStates; i++) {
				for (int j=0; j<numStates; j++) {
					double sum = 0.0;
					for(int s=0; s<numStates; s++)
						sum += (*ptr++) * eigValExp[s];
					branchMatrix(i, j) = (sum < 0.0) ? 0.0 : sum;
				}
			}

			returnMatrices.push_back(branchMatrix);
		}
	}
	else {
		for(int i = 0; i < branches.size(); i++){
			double v = branches[i];

			std::vector<std::complex<double>> ceigValExp;

			Matrix<double> branchMatrix(numStates, numStates);

			for (int s=0; s<numStates; s++)
				ceigValExp.push_back(exp(localComplexRateEigen.ceigenvalue[s] * v));

			const std::complex<double>* ptr = localComplexRateEigen.cc_ijk;
			for (int i=0; i<numStates; i++){
				for (int j=0; j<numStates; j++) {
					std::complex<double> sum = std::complex<double>(0.0, 0.0);
					for(int s=0; s<numStates; s++)
						sum += (*ptr++) * ceigValExp[s];
					branchMatrix(i, j) = (sum.real() < 0.0) ? 0.0 : sum.real();
				}
			}
		}
	}
}

// Be sure you want to delete!!
void TransitionProbability::deleteQ(const int index) {
	isComplex.erase(isComplex.begin() + index);
	rateEigen.erase(rateEigen.begin() + index);
	complexRateEigen.erase(complexRateEigen.begin() + index);
}

void TransitionProbability::deleteNQ(const int count) {
	for(int i = 0; i < count; i++){
		isComplex.pop_back();
		rateEigen.pop_back();
		complexRateEigen.pop_back();
	}
}