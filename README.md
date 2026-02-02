# About

This software samples from the conditional distribution of evolutionary selection (*dN/dS*) on protein sequences or taxa, given a multiple sequence alignment (MSA, in Nexus format) and a fixed tree topology. This conditional distribution (known as the posterior distribution) does not have a closed form, so we approximate it using Markov Chain Monte Carlo [MCMC](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling. The primary goal is to estimate how selection varies both across taxa and among sites within a protein. To achieve this, the model combines [Reversible-Jump MCMC](https://en.wikipedia.org/wiki/Reversible-jump_Markov_chain_Monte_Carlo) to infer heterogeneity in evolutionary regimes over time with a [Chinese Restaurant Process](https://en.wikipedia.org/wiki/Reversible-jump_Markov_chain_Monte_Carlo) prior to infer heterogeneity in regime evolution across sites. Something important to note about our software is that we do not support ambiguous codons currently - all incomplete codons will be treated as gaps (completely ambiguous data). This may not be desirable if your alignments have a lot of incomplete codons (e.g., ATN, GNT).

# Installation
...

# Usage
...

# Acknowledgements
Some source files (mainly within `src/core/`) originated from a codebase I developed with John Huelsenbeck, particularly those involving Eigen decomposition and linear algebra routines. As this project evolved beyond the original scope, these components were adapted and extended, though the modifications remain relatively minor.