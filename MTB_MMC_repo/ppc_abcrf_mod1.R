#' Posterior predictive checks for the best fitting genealogy model
#' Using scaled mutation rate theta equaling the generalized Watterson estimator
#' This just covers the simulations under the median posterior parameter
#' The input file has to be specified in line 77
#' It also allows input files with the same header but multiple data sets
#' i.e. gluing abc output together (with the header just once on top).

#' R packages
library(abcrf)
library(parallel)
library(psych)
library(gap)
library(phyclust)

mc1 <- 7 #cores used for parallel computation
nsim_ppc <- 10000 #Number of sims for the PPC

#' Observed data
load("TB_data_ABC_data_fasta.RData")

#' Seed (if wanted, used differently in the manuscript)
load("seeds_TB.RData")
seed1 <- seeds1[23*4+1]  
set.seed(seed1)



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



#' Diversity stats for PPC
divfun_check <- function(seq1){
  if (is.matrix(seq1)){
    out1 <- c(quant_hm_oc(seq1),mean_sd_oc(seq1),hammfun(seq1),
              f_nucdiv_S(spectrum01(seq1)),D_H(seq1)[1],allele_freqs(seq1)[-2])}
  else {out1 <- rep(NA,20)}
  names(out1) <-c("hm(O)",paste("O: qu",seq(.1,.9,.2)),"mean(O)","sd(O)",
                  paste("Ham: qu",seq(.1,.9,.2)),
                  "Nucl. div.","S","Taj.'s D",
                  paste("AF: qu",seq(.1,.9,.2))[-2])
  return(out1)}





#' Read in abc results (example: output of Rscript abcrf_tb.R 4 1)
#' THIS HAS TO BE ADJUSTED
abc_res <- read.table("restable_rev.txt",header = TRUE)


  

#' Best fitting model PPC 

#' Save PPC data and observed stats in the following objects
ppc_sims_fit1 <- vector("list",nrow(abc_res)) 
names(ppc_sims_fit1) <- abc_res$dataset
target1 <- ppc_sims_fit1

  
for (i in 1:nrow(abc_res)){ 
  if (is.na(abc_res$modelsel[i])){next}  
  i_rdata <- which(names(data_main)==abc_res$dataset[i])
  setwd("../general_scripts/")  
  prior1 <- prior_obs_s(data_main[[i_rdata]]$n_ind,
                        models=abc_res$modelsel[i],
                        nsimul=c(ifelse(abc_res$modelsel[i]==1,nsim_ppc,0),
                                 ifelse(abc_res$modelsel[i]==2,nsim_ppc,0),
                                 ifelse(abc_res$modelsel[i]==3,nsim_ppc,0),
                                 0,0,#just models 1,2,3,6
                                 ifelse(abc_res$modelsel[i]==6,nsim_ppc,0)),
                        ranges = list(rep(abc_res$postmed[i],2),
                                      rep(abc_res$postmed[i],2),
                                      rep(abc_res$postmed[i],2),
                                      0,0,0),# Kingman has always parameter 0, models 4,5 ignored 
                        s_obs = rep(data_main[[i_rdata]]$n_mut,2))
  
  setwd("../MTB_MMC_repo/")  
  
  clu1 <- makeForkCluster(nnodes = mc1)
  sims1 <- parApply(clu1,prior1,1,function(x){
    divfun_check(sim_seq(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],model = x[2]))})
  stopCluster(clu1)
  
  target1[[i]] <- as.data.frame(t(divfun_check(data_main[[i_rdata]]$seq_0_1)))
  ppc_sims_fit1[[i]] <- sims1
  
}


 

save(ppc_sims_fit1,target1,file="ppc_mod1_rev.RData")



