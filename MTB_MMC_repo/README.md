R Code for performing the ABC analyses from

"Multiple merger genealogies in outbreaks of Mycobacterium tuberculosis""  

F. Menardo, S. Gagneux and F. Freund

All scripts come with at least minimal documentation within, please check before running. Several scripts should be run via 'Rscript' on the command line, details are provided within the files. Settings etc are described in detail in our preprint. 

Subfolders

 * data_fasta/: Contains the FASTA files (only SNPs). Each files contains the SNPs for all strains and the reconstructed ancestral sequence of the dataset (named ">node*")
 
 * abc_res/: Empty folder, designated for output of abcrf_tb_popstruct.R, abcrf_tb_expd.R, abcrf_tb.R
 
 * resdata/: Empty folder, designated for output of multiple scripts

R scripts:

 * TB_datasets_for_ABC.RData: Contains SNP matrices for all data sets
 
 * abcrf_tb.R: Performs ABC model selection and parameter estimation for data sets from TB_datasets_for_ABC.RData. Run via command line (Rscript), details described within the file
   
 * abcrf_tb_popstruct.R: Performs ABC model selection and parameter estimation for Lee 2015 with population structure
   
 * abcrf_tb_expd.R: Performs ABC model selection and parameter estimation for Lee 2015 with exponential decline
   
 * serial_data.RData: Information about serial sampling for 3 data sets (Eldholm 2015, Lee 2015, Roetzer 2013)
 
 * serial_sampling_priorsim.R: Simulates serially sampled data under different coalescent models (Dirac, Beta, Kingman w. exp. growth)
   
 * sim_difft.R: Simulation wrapper for serial coalescents w. input from serial_data.RData
   
 * pe_serial_error.R, modsec_serial_error.R: Assesses errors in ABC parameter estimation/model selection of serial coalescents when random forest is built with standard (non-serial) coalescent simulations
   
 * seeds_TB.RData: Contains random seeds (only partially used)

* ppc_abcrf_mod1.R: Performs posterior predictive checks (based on posterior median coalescent parameter of the best fitting model) based on output from abcrf_tb.R

* ppc_plot.R: Plots posterior predictive checks based on output from ppc_abcrf_mod1.R 