# ABC model selection - full results 

## General information
Results are available as R objects

* Two folders: abc_res/ contains the ABC model selection results for different sets of statistics (output of abcrf_xfold_2rep.R in the parent folder distpaper_res/), while mismod_res/ contains the ABC results for sorting a set of simulations by using a random forest constructed to distinguish models from a different set of simulations (output of mA_abcrfvsmBC.R or setvsset_abcrf.R in the parent folder distpaper_res/).  


* Script oob.R: This script allows to directly extract the mean OOB error (OOB error for misclassifying one model class, averaged over all model classes) from the files listed below. In other words, for each replication, it shows the mean 
across the misclassification errors shown in each subfigure for our model selection error pictures. We used this script to get all the mean OOB error values in the manuscript. For usage, see the outcommented example in the script.
It also contains a function to extract the confusion matrix as a measure for model misclassification, averaged over both replications and 10 ABC runs for each replication.
 
## Contents of subfolder abc_res/



Objects that contain the loadings of the discriminant functions (not currently used)
 * lda_load_scen2_m12_n100_rep1.RData
 * lda_load_scen2_m12_n100_rep2.RData


All files starting with modsel_rep1, modsel_rep2 contain the results of 
repetitions 1 and 2 for the plots in the ms, so res_conf 


Structure of R object: contains 2 objects
 
 * res_conf: A list, each entry is again a list with the 10 replications of the ABC analysis with a specific set of statistics (sublists are named after the set of statistics used).
 * res_imp: A list of length 10. Each entry contains the variable importances for one ABC replication for a specific set of statistics (usually the full set of statistics)  

model numbers 1 = K+exp, 2 = Beta, 3 = Dirac, 4 = Beta-Xi, 5: Dirac+exp
              6= Kingman (7 not used here) 8= population structure 

Using AF+, r^2, O, Ham, Phy

 * modsel_rep*_m1235_af_n100.RData
 * modsel_rep1_m1235_af_n100_largeg.RData: As in the paper, but w. growth rates
                                           uniformly drawn from {1,10,50,100} 
                                           for Dirac+Exp
 * modsel_rep*_m123_af_n100.RData
 * modsel_rep*_m124_af_n100.RData
 * modsel_rep1_m12_all_n100af.RData
 * modsel_rep1_m12_af_TB.RData
 * modsel_rep1_scenario1.RData
 * modsel_rep1_m128_allmix1_n100af.RData: Mixed migration rates (see Appendix)
 * modsel_rep1_m128_allp5p5_n100af.RData: Sampling 50:50
 * modsel_rep1_m128_allp9p1_n100af.RData: Sampling 90:10
 * modsel_rep1_m13_af_n100.RData
 * modsel_rep1_m14_af_n100.RData
 * modsel_rep1_m16_af_n100.RData
 * modsel_rep1_m18_p5p5_allaf.RData
 * modsel_rep1_m18_p9p1_allaf.RData
 * modsel_rep1_m23_af_n100.RData
 * modsel_rep1_m24_af_n100.RData
 * modsel_rep1_m28_p5p5_allaf.RData
 * modsel_rep1_m28_p9p1_allaf.RData
 * modsel_rep1_m35_af_n100.RData
 * modsel_rep1_m12_all_n200af.RData
 * modsel_rep1_m12_all_n25af.RData


Using other stats
 * modsel_rep1_m128_p5p5_n100sfs.RData
 * modsel_rep1_m128_p9p1_n100sfs.RData
 * modsel_rep1_m12_all_n100fine.RData: More quantiles
 * modsel_rep1_m12_all_n100fold.RData
 * modsel_rep1_m12_all_n100lda.RData: With discriminant functions
 * modsel_rep1_m12_all_n100ohamsfs.RData 
 * modsel_rep1_m12_all_n100sfs.RData
 * modsel_rep1_m12_all_n200sfs.RData
 * modsel_rep1_m12_all_n25sfs.RData
 * modsel_rep1_m12_nos_n100.RData: No singletons
 * modsel_rep1_m12_rob_n100.RData: Scaled stats

## Contents of subfolder mismod_res/

 * modvsmod...
The files have been renamed to be more descriptive and now refer to the rows in Tables 8 and Table A2 in the manuscript. See the manuscript for more details.
 * simvssim_m123_af_n100_impc.RData: Contains a comparison between OOB errors and errors for misclassifying an independent test set for the comparison K+exp vs. Beta vs Dirac (output of setvsset_abcrf.R in the parent folder)

