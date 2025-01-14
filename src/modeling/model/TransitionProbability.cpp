#include <complex>
#include "TransitionProbability.hpp"
#include "core/Math.hpp"
#include <cstring>

TransitionProbability::TransitionProbability(const int nn)
    : numStates(122), numNodes(nn), probs1(), probs2() {

	/*
	probs[0] = new Matrix<double>*[2*numNodes*numCats];
    probs[1] = probs[0] + numNodes;

    for(int i = 0; i < numNodes*numCats; i++){
        probs[0][i] = new Matrix<double>(numStates, numStates, 0.0);
        probs[1][i] = new Matrix<double>(numStates, numStates, 0.0);
    }
	*/

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
	
	for(auto i : probs1){
		delete [] i;
	}
	for(auto i : probs2){
		delete [] i;
	}
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
	Matrix<double> P0 = (*this)(state, rate, node);
	if (!isComplex[rate])
		tiProbsEigens(v, P0, rateEigen[rate]);
	else
		tiProbsComplexEigens(v, P0, complexRateEigen[rate]);
}

/* This function calculates transition probabilities using
   complex eigenvalues and eigenvectors. */
void TransitionProbability::tiProbsComplexEigens(const double v, Matrix<double>& P, ComplexRateEigen& rE) {

	std::vector<std::complex<double>> ceigValExp;

	for (int s=0; s<numStates; s++)
		ceigValExp.push_back(exp(rE.ceigenvalue[s] * v));

	const std::complex<double>* ptr = rE.cc_ijk;
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
void TransitionProbability::tiProbsEigens(const double v, Matrix<double> &P, RateEigen& rE) {
	
	std::vector<double> eigValExp;

	for (int s=0; s<numStates; s++)
		eigValExp.push_back(exp(rE.eigenvalue[s] * v));

	double *ptr = rE.c_ijk;
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

			probs1.push_back(new Matrix<double>[numNodes]);
			probs2.push_back(new Matrix<double>[numNodes]);
			for(int j = 0; j < numNodes; j++){
				probs1.back()[j] = Matrix<double>(numStates, numStates, 0.0);
       	 		probs2.back()[j] = Matrix<double>(numStates, numStates, 0.0);
			}
		}
	}

	isComplex.shrink_to_fit();
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();

}

void TransitionProbability::updateQ(Matrix<double> Q, const int index) {
	isComplex[index] = eigens->update(Q, rateEigen[index], complexRateEigen[index]);
}

// Be sure you want to delete!!
void TransitionProbability::deleteQ(const int index) {
	isComplex.erase(isComplex.begin() + index);
	rateEigen.erase(rateEigen.begin() + index);
	complexRateEigen.erase(complexRateEigen.begin() + index);

	auto prob_it1 = probs1.begin() + index;
	delete [] *prob_it1;
	probs1.erase(prob_it1);

	auto prob_it2 = probs2.begin() + index;
	delete [] *prob_it2;
	probs2.erase(prob_it2);

	isComplex.shrink_to_fit();
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}

void TransitionProbability::deleteNQ(const int count) {
	for(int i = 0; i < count; i++){
		isComplex.pop_back();
		rateEigen.pop_back();
		complexRateEigen.pop_back();

		auto probs_it1 = std::prev(probs1.end());
		delete [] *probs_it1;
		probs1.pop_back();

		auto probs_it2 = std::prev(probs2.end());
		delete [] *probs_it2;
		probs2.pop_back();
	}
	
	isComplex.shrink_to_fit();
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}