R Code for performing the ABC analyses from

"Multiple merger genealogies in outbreaks of Mycobacterium tuberculosis""  

F. Menardo, S. Gagneux and F. Freund

All scripts come with at least minimal documentation within, please check before running. Several scripts should be run via 'Rscript' on the command line, details are provided within the files. Settings etc are described in detail in our preprint. 

Subfolders

 * data_call* : Contains the FASTA files (only SNPs) for two SNP calling procedures (minimum proportion % of reads have to show the variant). Each files contains the SNPs for all strains and the reconstructed ancestral sequence of the dataset (named ">node*")
 
 * all folder res*, folder sims: Empty folders in which our scripts write output 

 * plot_scripts_archive: Scripts for all ABC-related plots from the manuscript. Needs analysis to be rerun (or at least the extraction of the necessary results from the full ABC results). This is for documentation/archiving, scripts will not run directly.
 
 * BSP: Includes all raw data for Bayesian Skyline Plots from SNP matrices generated under Beta coalescents and settings close to data set Eldholm 2015 (simulated via sim_eldholm_clone.R), and one sample input file for Beast.

 * abc_results_preprint: Includes ABC results from all analyses. Folder names describe how analyses within were performed. Load R package abcrf for convenience when looking at the results. Also includes the results ppc_analysis1.RData for the posterior predictive checks. All objects within all subfolders but mmcg and serial contain the result of the model selection (pred_ms) and the parameter estimates under the Dirac coalescent, Beta coalescent and Kingman's n-coalescent with exponential growth for the specified data set.    
    * unig: Uniform prior on exponential growth
    * logg: Log uniform prior on exponential growth
    * beta0: Uniform prior on [0,2] for Beta coalescent
    * beta1: Uniform prior on [1,2] for Beta coalescent
    * largeth: Using larger prior range for scaled mutation rate theta
    * call* : Which data set is used (call75 or call90)
    * bsz: If an atom on alpha=1 is included for Beta coalescent (to include simulations of the Bolthausen-Sznitman coalescent)
    * mmcg: ABC for multiple merger coalescents with growth. Each file includes the confusion matrix (conf.matrix) and the model selection result (data2model). 
    * serial: ABC analysis under serial sampling (on simulated data), contains the results of the scripts pe_serial_error_log.R, modsec_serial_error_log.R.

R scripts/RData objects:

 * mtb_data_call75.RData, mtb_data_call90.RData: Contains SNP matrices for all data sets for both SNP calling settings
 
 * sim_all.R, abcrf_all.R: Performs various coalescent simulations (first script) and ABC analyses (model selection, parameter estimation) using these (second script) for the data sets in the files listed above. See the documentation within the two scripts.
   
 * abcrf_tb_popstruct.R: Performs coalescent simulation and ABC model selection and parameter estimation for Lee 2015 with population structure
   
 * abcrf_tb_expd.R: Performs coalescent simulation and ABC model selection and parameter estimation for Lee 2015 with exponential decline
 
 * abcrf_tb_mmcexp_savesims.R: Performs coalescent simulation and ABC model selection with added model classes Dirac and Beta coalescents in exponentially growing populations (and saves the simulations as RData object)
 
 * abcrf_tb_cont_largethetavar.R:  Performs coalescent simulation and ABC model selection and parameter estimation for the main ABC analysis, but with larger uncertaunty around the scaled mutation rate theta
 
 * serial_data.RData: Information about serial sampling for 3 data sets (Eldholm 2015, Lee 2015, Roetzer 2013)
 
 * serial_sampling_priorsim.R: Simulates serially sampled data under different coalescent models (Dirac, Beta, Kingman w. exp. growth)
   
 * sim_difft.R: Simulation wrapper for serial coalescents w. input from serial_data.RData
   
 * pe_serial_error_log.R, modsec_serial_error_log.R: Assesses errors in ABC parameter estimation/model selection of serial coalescents when random forest is built with standard (non-serial) coalescent simulations
   
 * seeds_TB.RData: Contains random seeds (only partially used)

* ppc_abcrf_mod1.R: Performs posterior predictive checks (based on posterior median coalescent parameter of the best fitting model, includes simulation) based on output from abcrf_tb.R

* prepare_abc_input.R: Produces mtb_data_call75.RData, mtb_data_call90.RData from the raw data.

* sim_eldholm_clone.R: Generates SNP matrices under Beta coalescents and settings close to data set Eldholm 2015