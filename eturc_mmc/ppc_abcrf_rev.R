#' Posterior predictive checks for the best fitting genealogy model
#' Using scaled mutation rate theta equaling the generalized Watterson estimator
#' This just covers the simulations under the median posterior parameter
#' The input file has to be specified in line 77
#' It also allows input files with the same header but multiple data sets
#' i.e. gluing abc output together (with the header just once on top).

#' R packages
#library(abcrf)
library(parallel)
library(psych)
library(gap)
library(phyclust)

mc1 <- 3#10 #cores used for parallel computation
nsim_ppc <- 75000#10000 #Number of sims for the PPC

#' ABC results
#' We use our weighted mean estimate for KM-exp and Beta
best_model <- c("red"=1,"green"=1,"lightblue"=1,"pink"=2,"kenya"=1)
load("resdata/abc_wmean_g_rev.RData")
best_param <- wmed_g
#best_param["pink"] <- 1.71 #As estimated via analyse_abc_2m_eturc.R

#' Necessary R code
source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/ext_fun_CREB.R")
source("../general_scripts/divstats.R")
source("../general_scripts/ext_fun_TajD_FWH.R")
source("../general_scripts/lambdacoal_sim.R")
source("../distpaper_res/construct_prior.R")

#' Further function wrappers


#' Wrappers for sequence simulation and diversity statistics
sim_seq <- function(nsamp1,theta1,coal_param=0,model=1){
  KM <- (model==1 & coal_param==0) | (model==2 & length(coal_param)==1 & coal_param==2) | (model==6)
  if (KM){raw_seq <- ms(nsam=nsamp1,opts=paste("-t",theta1))
  seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])}
  if (model==1 & !KM){
    raw_seq <- ms(nsam=nsamp1,opts=paste("-t",theta1,"-G",coal_param))
    seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])
  }
  if (model==2 & !KM){
    if (length(coal_param)==1){coal_param <- c(2-coal_param[1],coal_param[1])}
    if (all(coal_param == c(1,1))){seq1 <- bsz_seq_sim(nsamp1,theta1)} else {
      seq1 <- beta_seq_sim(nsamp1,coal_param[1],coal_param[2],theta1)}
  }
  if (model==3){
    seq1 <-  dirac_seq_sim(nsamp1,coal_param,theta1)
  }
  if (sum(seq1)==0){return(NA)}
  return(seq1)
}






  

#' Best fitting model PPC 

#' Save PPC data and observed stats in the following objects
#' Save PPC data and observed stats in the following objects
ppc_sims_fit1 <- list(vector("list",5),vector("list",5),vector("list",5),
                      vector("list",5),vector("list",5))
names(ppc_sims_fit1) <-  c("red","green","lightblue","pink","kenya")

target1 <- ppc_sims_fit1

  
for (lineage in c("green","red","lightblue","kenya")){ 
 load(paste0("input_data/ABCinput_",lineage,"_18scaff.RData"))
  for (scaff in 1:5){
  setwd("../general_scripts/")  
  prior1 <- prior_obs_s(data_scaffs[[scaff]]$n_ind,
                        models=best_model[lineage],
                        nsimul=c(ifelse(best_model[lineage]==1,nsim_ppc,0),
                                 ifelse(best_model[lineage]==2,nsim_ppc,0),
                                 ifelse(best_model[lineage]==3,nsim_ppc,0),
                                 0,0,#just models 1,2,3,6
                                 ifelse(best_model[lineage]==6,nsim_ppc,0)),
                        ranges = list(rep(best_param[lineage],2),
                                      rep(best_param[lineage],2),
                                      rep(best_param[lineage],2),
                                      0,0,0),# Kingman has always parameter 0, models 4,5 ignored 
                        s_obs = rep(data_scaffs[[scaff]]$n_mut,2))
  
  setwd("../eturc_mmc/")  
  
  if (data_scaffs[[scaff]]$n_ind<30){lump1 <- FALSE} else {lump1 <- TRUE} 
  divfun_fold_check <- function(seq1,n_ind){
    if (lump1){l1 <- 17} else {l1 <- 2 + ceiling((n_ind-1)/2)}
    if (is.matrix(seq1)){
      out1 <- c(f_nucdiv_S(spectrum01(seq1))[1],
                fnl_sfs(spectrum01(seq1),TRUE,FALSE,lump1,15),
                D_H(seq1)[1])} else {
                  out1 <- rep(NA,l1)}
    if (lump1){
    names(out1) <-c("Nucleotide diversity",
                    paste0("fS",1:14),"fS15+","Tajima's D")} else {
            names(out1) <-c("Nucleotide diversity",
                             paste0("fS",1:(l1-2)),"Tajima's D")      
                    }
    return(out1)} 
  clu1 <- makeForkCluster(nnodes = mc1)
  sims1 <- parApply(clu1,prior1,1,function(x){
    divfun_fold_check(n_ind = x[1], seq1 = 
    sim_seq(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],model = x[2]))})
  stopCluster(clu1)
  
  target1[[lineage]][[scaff]] <- as.data.frame(t(divfun_fold_check(data_scaffs[[scaff]]$seq_0_1,
                                                           data_scaffs[[scaff]]$n_ind)))
  ppc_sims_fit1[[lineage]][[scaff]] <- sims1
  
}}


 

save(ppc_sims_fit1,target1,file="ppc_eturc_TajD_rev2.RData")



