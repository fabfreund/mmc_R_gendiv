# ABC model selection - full results 

Results are available as R objects

Objects that contain the loadings of the discriminant functions (not currently used)
 * lda_load_scen2_m12_n100_rep1.RData
 * lda_load_scen2_m12_n100_rep2.RData


oob.R (convenience copy from oob.R in parent folder)
 * This script allows to directly extract the mean OOB error (OOB error for misclassifying 
  one model class, averaged over all model classes) from the files listed below. In other words, for each replication, it shows the mean 
across the misclassification errors shown in each subfigure for our model selection error
pictures. We used this script to get all the mean OOB error values in the manuscript.
For usage, see the outcommented example in the script.

----------------------------------

All files starting with modsel_rep1, modsel_rep2 contain the results of 
repetitions 1 and 2 for the plots in the ms, so res_conf 


Structure of R object: contains 2 objects
 
 * res_conf: A list, each entry is again a list with the 10 replications of the ABC analysis with a specific set of statistics (sublists are named after the set of statistics used).
 * res_imp: A list of length 10. Each entry contains the variable importances for one ABC replication for a specific set of statistics (usually the full set of statistics)  

model numbers 1 = K+exp, 2 = Beta, 3 = Dirac, 4 = Beta-Xi, 5: Dirac+exp
              6= Kingman (7 not used here) 8= population structure 

Using AF+, r^2, O, Ham, Phy

 * modsel_rep*_m1235_af_n100.RData
 * modsel_rep*_m123_af_n100.RData
 * modsel_rep*_m124_af_n100.RData
 * modsel_rep1_m12_all_n100af.RData
 * modsel_rep1_m12_af_TB.RData
 * modsel_rep1_scenario1.RData
 * modsel_rep1_m128_mix1_n100af.RData: Mixed migration rates (see Appendix)
 * modsel_rep1_m128_p5p5_n100af.RData: Sampling 50:50
 * modsel_rep1_m128_p9p1_n100af.RData: Sampling 90:10
 * modsel_rep1_m13_af_n100.RData
 * modsel_rep1_m14_af_n100.RData
 * modsel_rep1_m16_af_n100.RData
 * modsel_rep1_m18_p5p5_n100af.RData
 * modsel_rep1_m18_p9p1_n100af.RData
 * modsel_rep1_m23_af_n100.RData
 * modsel_rep1_m24_af_n100.RData
 * modsel_rep1_m28_p5p5_n100af.RData
 * modsel_rep1_m28_p9p1_n100af.RData
 * modsel_rep1_m35_af_n100.RData
 * modsel_rep1_m12_all_n200af.RData
 * modsel_rep1_m12_all_n25af.RData
 * modsel_rep1_m12_smear2_lowg_n100.RData: 
 * modsel_rep1_m12_wattest_n100.RData:


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
