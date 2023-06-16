This folder contains the input data, analysis scripts and output files for the ABC analysis from Vidal-Villarejo et al. "Population history of the Northern corn leaf blight fungal pathogen *Setosphaeria turcica* in Europe" (2020, preprint on biorxiv). Essentially just an application of the methodology used for *Mycobacterium tuberculosis* (see other folder in this repository), but without polarized SNP data. 

Subfolders

  * Folder input_data/: Contains the SNP data of the 18 largest scaffolds for the five clusters Big Clonal (in file names =red), Small Clonal (=green), French Clonal (=lightblue), Diverse (=pink) and Kenyan (=kenya)
  * Folder res/: Contains the results of model selection and parameter estimation between/within Kingman's n-coalescent with exponential growth, Beta n-coalescents and Dirac n-coalescents. Output of script abcrf_eturc_m123.R.
  Folder res2m/: Contains the results of model selection and parameter estimation between/within Kingman's n-coalescent with exponential growth and the better fitting model from Beta n-coalescents and Dirac n-coalescents (based on the three-fold ABC model selection above). Output of script abcrf_eturc_2m.R.
  * Folder sims/: Empty folder, will contain the simulations used to construct the random forests for the ABC for model selection between Kingman's n-coalescent, Beta n-coalescents and Dirac n-coalescents with exponential growth (for each model) if script abcrf_eturc_mmcexp.R is run. Simulations not included due to size, but reproducible. 
  * Folder resdata/: Contains the ABC results for model selection between Kingman's n-coalescent, Beta n-coalescents and Dirac n-coalescents with exponential growth (for each model). Output of script abcrf_eturc_mmcexp.R.

Files

 * The ABC and simulations were done by running the scripts abcrf_eturc_m123.R, abcrf_eturc_2m.R and abcrf_eturc_mmcexp.R via Rscript (details are provided as comments within each script) 
 * ppc_abcrf.R simulates genetic diversity statistics under the posterior median growth rate fitted to each cluster and also provides the actually observed values (output is ppc_eturc_TajD_rev2), ppc_plot.R then performs the graphical posterior predictive checks (output ppc_eturc.pdf).
 * extract_mmcexpres.R: Convenience scripts to extract key results from the result files for the MMC-exp model selection (not really interesting tbh).
 
 