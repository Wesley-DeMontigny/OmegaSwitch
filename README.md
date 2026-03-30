# OmegaSwitch

This software samples from the conditional distribution of evolutionary selection (*dN/dS*) on protein sequences or taxa, given a multiple sequence alignment (MSA, in Nexus format) and a fixed tree topology. This conditional distribution (known as the posterior distribution) does not have a closed form, so we approximate it using Markov Chain Monte Carlo [MCMC](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling. The primary goal is to estimate how selection varies both across taxa and among sites in a sequence alignment. To achieve this, the model combines [Reversible-Jump MCMC](https://en.wikipedia.org/wiki/Reversible-jump_Markov_chain_Monte_Carlo) to infer heterogeneity in evolutionary regimes over time with a [Dirichlet process](https://en.wikipedia.org/wiki/Dirichlet_process) prior to infer heterogeneity in regime evolution across protein residues. We utilize ancestral sampling over the course of the MCMC algorithm to produce posterior distributions of *dN/dS* at each entry in the alignment. For example, this is the posterior mean *dN/dS* acting on a particular alpha-helix in alpha-hemoglobin (grey indicates a low *dN/dS*, while red indicates a high one.).

![Alpha Globin](/publication_analyses/globin_analysis/dpcmm/alpha_solv_dpcmm.png)

This software implements three different models: M0, CMM, and DPCMM. The M0 model is a typical codon phylogenetic model. The CMM model is a Codon Markov Modulated model, which uses RJ-MCMC to select the number of hidden evolutionary regimes that Markov-Modulated CTMC can swap between. The DPCMM model takes this further, using RJ-MCMC to infer the number of hidden evolutionary regimes and modeling the data as an infinite mixture of these processes trhough the dirichlet process. This final model is an incredibly computationally intensive model, and to ensure proper sampling from the posterior, we recommend using our Metropolis-Coupled MCMC (MCMCMC) sampler. In our experience, a beta of about 0.8 results in around a 10% acceptance rate for MCMCMC chain swaps under this model. Additionally, inferring the mixture over Markov-modulated processes, while simultaneously inferring the number of regimes can result in some non-identifiability between the two parameters. However, in this case the model should still recover correct posterior distributions for the *dN/dS* across the tree.

If you want to disable RJ-MCMC over the hidden regime count, you can provide `-fixedRegimes` during `CMM` or `DPCMM` inference. This initializes the model at the requested number of regimes and prevents any regime-count RJ proposals for the rest of the run. Valid bounds are `1-5` for `CMM` and `1-3` for `DPCMM`.

Something important to note about our software is that we do not support ambiguous codons currently - all incomplete codons will be treated as gaps (completely ambiguous data). This may not be desirable if your alignments have a lot of incomplete codons (e.g., ATN, GNT).



# Installation
...

# Usage

```
Inference Input/Output:
   * -nexus             : Input nexus file containing the nculeotide alignment.
   * -treeOut           : The output file name for the tree trace.
   * -mcmcOut           : The output file name for the bulk of the MCMC trace, excluding the tree and DPP parameters.
   * -dppOut            : The output file name for the DPP parameters.
   * -tipsOut           : The output file name for the reconstructed tip dNdS ratios.
   * -ancestralStatesOut: The output file name for the all ancestral dNdS ratios.
   * -tree              : The NEWICK string corresponding to the fixed tree you wish to analyze.
   * -simulationOutput  : The output file name for the true simulation parameters.
   * -threads           : The number of threads to use during the analysis.

Inference Model and Simulation:
   * NOTE: By default, this software will run the M0 model.
   * -M0                : Do inference under a normal codon phylogenetic model.
   * -CMM               : Do inference under the reversible jump Markov-modulated model.
   * -DPCMM             : Do inference with the reversible-jump DPP model.
   * -simulateM0        : Directs the program to simulate under M0 and test against the selected inference model.
   * -simulateCMM       : Directs the program to simulate under the CMM model and test against the selected inference model.
   * -simulateDPCMM     : Directs the program to simulate under the CMM-DPP model and test against the selected inference model.
   * -numSimulations    : The number of simulations to do inference under.

Model Parameters:
   * -treeMean          : Mean for the tree length prior.
   * -treeSD            : SD for the tree length prior.
   * -omegaLambda       : Rate parameter for the nonsynonymous mutation rate's exponential prior.
   * -kLambda           : Rate parameter for the transition/transversion rate's exponential prior.
   * -rLambda           : Rate parameter for the matrix-swapping rate's exponential prior.
   * -expectedCat       : The number of expected categories for the DPP.
   * -fixedRegimes      : Fix the number of hidden regimes for CMM or DPCMM and disable RJ-MCMC over regime count.
   * -mcmcmcBeta        : Inverse temperature for the auxiliary heated chain (setting this enables MCMCMC).

Sampling Options:
   * -numIter           : The number of iterations for the MCMC.
   * -printFreq         : How often to output the MCMC state to the screen.
   * -sampleFreq        : How often to ouput the MCMC state to log files.
   * -burnInIter        : The number of iterations for the burn-in.
   * -tuneFreq          : How often to tune the acceptance rate of the MCMC moves during the burn-in.
   * -bayesOpt          : (Experimental) The number of iterations to run Bayesian optimization on the MCMC moves after the burn-in.
   * -bayesOptFreq      : (Experimental) How often to sample the trace for Bayesian optimization.
   * -numGibbs          : How many Gibbs updates to perform on the DPP partitions.
   * -mcmcmcSwapFreq    : How many MCMC iterations between attempted swaps in MCMCMC.
   * -treeWeight        : How often to propose a move on the tree.
   * -kWeight           : How often to propose a move on the K parameter.
   * -rWeight           : How often to propose a move on the R parameter.
   * -rjWeight          : How often to propose a move on the dimension of the Markov modulated model.
   * -stationaryWeight  : How often to propose a move on the stationary distribution.
   * -dppWeight         : How often to propose a move on the DPP partitions.
   * -omegaWeight       : How often to propose a move on the omega parameters.
```

# Citing Us
...

# Acknowledgements
Some source files (mainly within `src/misc/`) originated from a codebase I worked on with Dr. John Huelsenbeck, particularly those involving Eigen decomposition and linear algebra routines. As this project evolved beyond the original scope, these components were adapted and extended, though the modifications remain relatively minor.
