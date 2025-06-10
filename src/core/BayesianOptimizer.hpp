#ifndef BAYESIAN_OPTIMIZER_HPP
#define BAYESIAN_OPTIMIZER_HPP

#include <vector>

class BayesianOptimizer {
    public:
        BayesianOptimizer(void)=delete;
        BayesianOptimizer(int nP, int s);
        ~BayesianOptimizer();

        static double objective(std::vector<double> r);
        std::vector<double> maximizeAcquisition(std::vector<double> s, double o); // TODO: Self explainatory
        std::vector<double> getMaximum(); // TODO: Final thing to implement
        void updateGaussianProcess(); // TODO: Implement Cholesky decomposition instead of LU. Update choleskyFactor and alpha. Update hyperparameters
        void initializeHyperparameters(); // TODO: Use latin hypercubes to set initial hyperparameters?
        void updateHyperparameters(); // TODO: Maximize the likelihood of hyperparameters

        void registerSample(std::vector<double> s) {samples.push_back(s);}
        void setScale(std::vector<double> s) {scale = s;}
    private:
        int numParams;
        int iterationsPerSample;

        std::vector<std::vector<double>> samples; // The MCMC parameters associated with a particular batch of MCMC iterations
        std::vector<double> objectives; // The values of the objective function for each of the samples

        std::vector<double> hyperparams; // The length scale of the ARD kernel
        std::vector<double> scale; // Gives us a good judge of how to box off our search.

        Matrix<double> choleskyFactor;
        std::vector<double> alpha;

        static double autocorrelationScore(const std::vector<double>& r, int end); //Mahendran et al. 2010
        double kernel(std::vector<double> a, std::vector<double> b); // ARD kernel
};

#endif