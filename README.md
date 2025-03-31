# strepto_demography

This repository contains code to run the ABC parameter estimation for samples from *S. pneumoniae* and *S. mitis* as well as the subsequent Ne estimation approach.
This is linked to the article "Long-term evolution of Streptococcus mitis and Streptococcus pneumoniae leads to higher genetic diversity within rather than between human populations" from Davison et al. (doi: https://doi.org/10.1371/journal.pgen.1011317).

Please see the infile documentation for each script for details how to use it. 
 
 * The two scripts expg_estim_... provide growth parameter estimation via Approximate Bayesian Computation (ABC). The observed data is provided by the RData objects (containing observed diversity statistics) in the subfolder /data. To provide training data for the ABC, we use a simulation approach that samples SNPs from expected genome-wide SNP allele frequency spectra (site frequency spectra, SFS). While the expected SFS can be computed (details described in the infile documentation in the R scripts), they are provided in /data for the analysis from our article).
 * The script Ne_comp_smitispneu.R estimates the effective population size N_e as described in the published article. For this, you need to compile CREBiepg.c which is called by ext_fun_CREB.R (details documented in ext_fun_CREB.R).
