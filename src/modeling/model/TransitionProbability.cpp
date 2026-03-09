#include <complex>
#include "TransitionProbability.hpp"
#include "misc/Math.hpp"
#include <cstring>

TransitionProbability::TransitionProbability(const int nn, const int ss)
    : numStates(ss), numNodes(nn) {

	Matrix<double> Q(numStates, numStates, 0.0);

	eigens = new EigenSystem;
	
	allocateQ(1);
	updateQ(Q, 0);
	accept();
}

TransitionProbability::~TransitionProbability(void) {
	delete eigens;
	
	for(auto i : probs[0]){
		delete [] i;
	}
	for(auto i : probs[1]){
		delete [] i;
	}
}


void TransitionProbability::accept(void) {
	isComplex[1] = isComplex[0];
	for(int i = 0; i < isComplex[0].size(); i++){
		if(!isComplex[0][i]){
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
	isComplex[0] = isComplex[1];
	for(int i = 0; i < isComplex[1].size(); i++){
		if(!isComplex[0][i]){
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
	if (!isComplex[0][rate])
		computeProbs(v, P0, rateEigen[rate]);
	else
		computeProbs(v, P0, complexRateEigen[rate]);
}

void TransitionProbability::setProbs(const int state, const int rate, const int node, const int stateSubset, const double v) {
	Matrix<double> P0 = (*this)(state, rate, node);
	if (!isComplex[0][rate])
		computeProbs(v, P0, rateEigen[rate], stateSubset);
	else
		computeProbs(v, P0, complexRateEigen[rate], stateSubset);
}

void TransitionProbability::computeProbs(const double v, Matrix<double>& P, ComplexRateEigen& rE) {

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


void TransitionProbability::computeProbs(const double v, Matrix<double> &P, RateEigen& rE) {
	
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

void TransitionProbability::computeProbs(const double v, Matrix<double>& P, ComplexRateEigen& rE, const int stateSubset) {

	std::vector<std::complex<double>> ceigValExp;

	for (int s=0; s<stateSubset; s++)
		ceigValExp.push_back(exp(rE.ceigenvalue[s] * v));

	const std::complex<double>* ptr = rE.cc_ijk;
	for (int i=0; i<stateSubset; i++)
		{
		for (int j=0; j<stateSubset; j++) 
			{
			std::complex<double> sum = std::complex<double>(0.0, 0.0);
			for(int s=0; s<stateSubset; s++)
				sum += (*ptr++) * ceigValExp[s];
			P(i, j) = (sum.real() < 0.0) ? 0.0 : sum.real();
			}
		}
}

void TransitionProbability::computeProbs(const double v, Matrix<double> &P, RateEigen& rE, const int stateSubset) {
	
	std::vector<double> eigValExp;

	for (int s=0; s<stateSubset; s++)
		eigValExp.push_back(exp(rE.eigenvalue[s] * v));

	double *ptr = rE.c_ijk;
	for (int i=0; i<stateSubset; i++) 
		{
		for (int j=0; j<stateSubset; j++) 
			{
			double sum = 0.0;
			for(int s=0; s<stateSubset; s++)
				sum += (*ptr++) * eigValExp[s];
			P(i, j) = (sum < 0.0) ? 0.0 : sum;
			}
		}
}

void TransitionProbability::allocateQ(int size){
	if(size > isComplex[0].size()) {
		for(int i = 0, num = size - isComplex[0].size(); i < num; i++){
			isComplex[0].push_back(false);
			rateEigen.push_back(RateEigen(numStates));
			complexRateEigen.push_back(ComplexRateEigen(numStates));

			probs[0].push_back(new Matrix<double>[numNodes]);
			probs[1].push_back(new Matrix<double>[numNodes]);
			for(int j = 0; j < numNodes; j++){
				probs[0].back()[j] = Matrix<double>(numStates, numStates, 0.0);
       	 		probs[1].back()[j] = Matrix<double>(numStates, numStates, 0.0);
			}
		}
	}
}

void TransitionProbability::updateQ(Matrix<double> Q, const int index) {
	isComplex[0][index] = eigens->update(Q, rateEigen[index], complexRateEigen[index]);
}


void TransitionProbability::deleteQ(const int index) {
	isComplex[0].erase(isComplex[0].begin() + index);
	rateEigen.erase(rateEigen.begin() + index);
	complexRateEigen.erase(complexRateEigen.begin() + index);

	auto prob_it1 = probs[0].begin() + index;
	delete [] *prob_it1;
	probs[0].erase(prob_it1);

	auto prob_it2 = probs[1].begin() + index;
	delete [] *prob_it2;
	probs[1].erase(prob_it2);
}

void TransitionProbability::deleteNQ(const int count) {
	for(int i = 0; i < count; i++){
		isComplex[0].pop_back();
		rateEigen.pop_back();
		complexRateEigen.pop_back();

		auto probs_it1 = std::prev(probs[0].end());
		delete [] *probs_it1;
		probs[0].pop_back();

		auto probs_it2 = std::prev(probs[1].end());
		delete [] *probs_it2;
		probs[1].pop_back();
	}
}