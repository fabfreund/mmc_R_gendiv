# R code to replicate results in 'Distinguishing coalescent models - which statistics matter most?'

IN "BETA TESTING", comments very welcome! 

The general approach is to first run one or multiple of the sim_.. scripts (in folder scripts_sim, sim_m8_all_n.R should be run via Rscript with a further input file, see sim_m8_all_n.R), then run abcrf_xfold_2rep.R with a parameter script. 
Simulation output (with one replication) is written into the folders sims_rep1/ and sims_rep2/. ABC results are written into folder abc_res/.  
Then, the ABC results can be plotted or summarised with additional scripts, see below.

Further documentation can also be found within the R scripts.


IMPORTANT: 

* Each simulation script comes with a specific set of diversity stats simulated, please adjust the sets of statistics to your needs. The specific setting from our paper simulated by each script should be recognizable from the file name.
* The ABC script abcrf_xfold_2rep.R performs 10 ABC runs on each replication for the simulated scenarios from our preprint (for each of several subsets of statistics). The specific analysis is specified in an input parameter script (which are in the subfolder scripts_params_abc). In order to perform the ABC, you need to first simulate as described above, and, where necessary, concatenate simulation runs. Concatenation is done by row-binding the prior matrices and column-binding the simulations within each replication and save each new object as <FILENAME>.RData in sims_rep1/ resp. sims_rep2. This may involve some more fine-tuning as choosing a subset of statistics. Script gluesim_m123.R provides an example.
* If you need Watterson's estimator for Kingman's n-coalescent with exponential growth, you need to first comile the C code general_scripts/CREBiepg.c. See general_scripts/extfun.R for more information.

File description (within this main folder)

* construct_prior.R: Wrapper script to produce priors for model selection between classes of Lambda- and Xi-coalescents 
* divfunwrappers.R: Wrapper script to jointly compute different diversity statistics 
* abcrf_xfold_2rep.R: General ABC analysis script (using R package abcrf). To be run with R script, details are comments within the script. 
* plot_abc_paper.R: Plots as in our preprint from output of abcrf_xfold_2rep.R. **Needs to be adjusted by hand to match your output** 
* oob.R: Summarises OOB errors of ABC analyses with abcrf_xfold_2rep.R. Currently simultaneously for two ABC analyses (or twice the same...)
* mA_abcrfvsmBC.R: Takes simulations from at least 3 model classes as input. Takes simulations from one model class and performs an ABC model selection for these simulations between the other models. 
* misclass_all_abc_m123.R: Comparison of different ABC methods on simulations of K+exp, Beta and Dirac coalescent classes. Test data needs to be simulated before, as described above (from all three models, to match our results run sim_m12_all_n100.R,  sim_m3_af_n100.R and then gluesim_m123.R) 
* gluesim_m123.R: Example of how we combined different simulation runs.
* folder raw_results_as_in_paper: Includes all ABC results reported in the paper as R objects. File name gives scenario, models involved and, if necessary, sample size. Coalescent models are abbreviated: 1 = Kingman w. exp growth, 2 = Beta, 3 = Dirac, 4 = BetaXi, 5 = Dirac+exp 6 = Kingman, 8 = Kingman coalescent with population structure and exponential growth. More details can be found in names2scenarios.txt within the folder
