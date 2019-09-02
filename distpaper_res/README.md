Currently: The folder contains sample code needed to replicate our analysis from 'Distinguishing coalescent models - which statistics matter most?'
THIS IS UNDER CONSTRUCTION. 

The general approach is to first run one or multiple of the sim_.. scripts, then run abcrf_xfold_2rep.R with a parameter script. 
Simulation output (with one replication) is written into the folders sims_rep1/ and sims_rep2/. ABC results are written into folder abc_res.  
Then, the ABC results can be plotted or summarised with additional scripts, see below.

Further documentation can also be found within the R scripts.


IMPORTANT: 
* Before running simulations and ABC model selections scripts, make new folders sims_rep1/, sims_rep2/ and/or abc_res within this folder
* Each simulation script comes with a specific set of diversity stats simulated, please adjust the sets of statistics to your needs.
* The ABC script abcrf_xfold_2rep.R performs 10 ABC runs on each replication for each of several subsets of statistics, which are specified in an input parameter script.
  Best have a look at one of the example scripts.
* If you need Watterson's estimator for Kingman's n-coalescent with exponential growth, you need to first comile the C code general_scripts/CREBiepg.c. See general_scripts/extfun.R for more information.

File description
* construct_prior.R: Wrapper script to produce priors for model selection between classes of Lambda- and Xi-coalescents 
* divfunwrappers.R: Wrapper script to jointly compute different diversity statistics 
* sim_m12_all_n*.R, sim_m12_all.R: Simulate (many) diversity statistics under Kingman's n-coalescent with exponential growth and under Beta-n-coalescents from prior distributions (Scenario 2 from the paper, n=*). 
* sim_m3_af_n100.R: Simulate diversity statistics under Dirac n-coalescents from prior distribution (Scenario 2, n=100)
* sim_m12_TB_n100.R: Simulate diversity statistics under Kingman's n-coalescent with exponential growth and under Beta-n-coalescents from prior distributions (Scenario 3 from the paper)
* sim_scenario1.R: Simulate diversity statistics in Scenario 1 (replication of Kato et al., Royal Society Open Science, 2017)
* abcrf_xfold_2rep.R: General ABC analysis script (using R package abcrf). To be run with R script and an input parameter script as the following
   * scenario2_sfsl.R: Using lumped SFS and further summary statistics
   * scenario2_fold.R: Using fSFS and further summary statistics that don't need polarized data
   * scenario3.R: ABC for Scenario 3 from the paper
   * kato_rep.R: ABC for Scenario 1
* plot_abc_paper.R: Plots as in our preprint from output of abcrf_xfold_2rep.R. **Needs to be adjusted by hand to match your output** 
* oob.R: Summarises OOB errors of ABC analyses with abcrf_xfold_2rep.R. Currently simultaneously for two ABC analyses (or twice the same...)
* mA_abcrfvsmBC.R: Takes simulations from at least 3 model classes as input. Takes simulations from one model class and performs an ABC model selection for these simulations between the other models. 
