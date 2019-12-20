#' ABC-RF on Mycobacterium tuberculosis data (sub)sets
#' Run via R script, two integer arguments 
#' Argument1: gives the list
#' position of the data set within the RData object TB_datasets_for_ABC.RData
#' Argument2: Replication run (changes seed and output name)


args <- commandArgs(TRUE)
i <- as.integer(args[1])



rep1 <- as.integer(args[2])

#' To make results reproducible, only partially applied
#' for the preprint, thus commented out
#load("seeds_TB.RData")
#seed1 <- seeds1[29*(rep1-1)+i]
#set.seed(seed1)

#' R packages 
library(abcrf)
library(parallel)
library(psych)
library(gap)
library(phyclust)

#' Scripts from this repo
source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/ext_fun_CREB.R")
source("../general_scripts/misc_scripts.R")
source("../general_scripts/lambdacoal_sim.R")
source("../general_scripts/divstats.R")
source("../distpaper_res/construct_prior.R")

#' Load the TB data
load("TB_datasets_for_ABC.RData")

#' number of cores for parallel computation
mc1 <- 7#16


#' Simulation parameter for model classes for the ABC
gran1 <- seq(0,5000,2) 
gran2 <- seq(0,20000,5) #same 
f1 <- function(x){switch(x,gran1,gran2)}
growthrange0 <- lapply(c(1,2,1,1,2,1,2,2,1,1,1,2,2,1,1,1,2,2,2,1,
                         1,1,1,1,2,2,2,2,2),f1)
betarange0 <- seq(1,1.975,0.025)
diracrange0 <- seq(0.025,0.975,0.025)
n0 <- 125000 # Number of simulations per model class

#' Build prior distribution in all model classes (including standard Kingman)
#' model 1: Kingman+exp, 2: Beta, 3: Dirac, 6: Kingman
setwd("../general_scripts/")
prior1 <- prior_obs_s(data_main[[i]]$n_ind,models=c(1,2,3,6),nsimul=c(n0,n0,n0,0,0,n0),
                      ranges = list(growthrange0[[i]],
                                    betarange0,diracrange0,NULL,NULL,0),
                      s_obs = rep(data_main[[i]]$n_mut,2))

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
setwd("../MTB_MMC_repo/")


#' Function to simulate under prior distribution (builds on coalescent
#' simulation, see ../general_scripts)
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

#' Function to compute genetic diversity statistics (summary statistics for ABC)
#' Option private: Kicks out private mutation (not used in the manuscript)
#' builds on other R scripts from the github repo
divfun_or2af <- function(seq1,private=TRUE){
  if (is.matrix(seq1)){
    log1 <- rep(TRUE,ncol(seq1))
    if (!private){log1 <- (colSums(seq1)>1)}
    if (sum(log1)==0){out1 <- rep(NA,20)} else {
      seq1 <- as.matrix(seq1[,log1])
      out1 <- c(quant_hm_oc(seq1),mean_sd_oc(seq1),
            r2fun(seq1),f_nucdiv_S(spectrum01(seq1)),allele_freqs(seq1))}}
  else {out1 <- rep(NA,20)}
        names(out1) <-c("o_hm",paste0("o_q",seq(1,9,2)),"o_mean","o_sd",
                        paste0("r2_q",seq(1,9,2)),"nucdiv","S",
                        paste0("AF",seq(1,9,2)))
  return(out1)}

#' Alternatively (large sample, many mutations), use different stats w/o r^2
#' This is not used when you run this script as is,
#' if you want to switch to these stats, change l. 127 and l. 143
divfun_few <- function(seq1,private=TRUE){
  if (is.matrix(seq1)){
    log1 <- rep(TRUE,ncol(seq1))
    if (!private){log1 <- (colSums(seq1)>1)}
    if (sum(log1)==0){out1 <- rep(NA,19)} else {
      seq1 <- as.matrix(seq1[,log1])
      out1 <- c(quant_hm_oc(seq1),mean_sd_oc(seq1),
                f_nucdiv_S(spectrum01(seq1)),allele_freqs(seq1,
                                                          quant_v = seq(.1,.9,.1)))}}
  else {out1 <- rep(NA,19)}
  names(out1) <-c("o_hm",paste0("o_q",seq(1,9,2)),"o_mean","o_sd",
                  "nucdiv","S",
                  paste0("AF_q",seq(1,9,1)))
  return(out1)}


#' Parallelized coalescent simulations

clu1 <- makeForkCluster(nnodes = mc1)
sims1 <- parApply(clu1,prior1,1,function(x){
  divfun_or2af(sim_seq(nsamp1 = x[1],theta1 = x[5],
                       coal_param = x[3],model = x[2]))})
stopCluster(clu1)

model1 <- prior1[,2]
coal_p <- prior1[,3]

#' For safety: Remove simulations with NA (due to overflow)
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})

sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
coal_p <- coal_p[!(bad_cols)]

sumstat1 <- t(sims1)
rownames(sumstat1) <- NULL
target1 <- as.data.frame(t(divfun_or2af(data_main[[i]]$seq_0_1)))
 

#' Parameter estimation Kingman+exp
param1 <- coal_p[model1==1]
sumstat2 <- data.frame(param1,sumstat1[model1==1,])
growth_rf <- regAbcrf(param1~.,
                      data=sumstat2,
                      paral = TRUE)
growthfit_rf <- predict(growth_rf,target1,
                        sumstat2,paral=TRUE,
                        quantiles = c(0,0.025,.05,.95,.975))

#' Parameter estimation Beta
param1 <- coal_p[model1==2]
sumstat2 <- data.frame(param1,sumstat1[model1==2,])  
beta_rf <- regAbcrf(param1~.,
                    data=sumstat2,
                    paral = TRUE)
betafit_rf <- predict(beta_rf,target1,sumstat2,paral=TRUE,
                      quantiles = c(0,0.025,.05,.95,.975))

#' Parameter estimation Dirac
param1 <- coal_p[model1==3]
sumstat2 <- data.frame(param1,sumstat1[model1==3,])   
dirac_rf <- regAbcrf(param1~.,
                     data=sumstat2,
                     paral = TRUE)
diracfit_rf <- predict(dirac_rf,target1,sumstat2,
                       paral=TRUE,quantiles = c(0,0.025,.05,.95,.975))



model1 <- as.factor(model1)
sumstat1a <- data.frame(model1,t(sims1)) 
#' Remove summary statistics with no variation within a model class
good_cols <- lapply(sumstat1a[,-1],
                    function(v){sum(tapply(v,sumstat1a[,1],function(x){range(x)[1]==range(x)[2]}))>=1}) #added that single no var classes already cause problems...
good_cols <- !(c(FALSE,unname(unlist(good_cols))))

#' ABC model selection
rf_ms <- abcrf(model1~.,data=sumstat1a[,good_cols],paral = TRUE)
pred_ms <- predict(rf_ms,target1[good_cols[-1]],sumstat1a[,good_cols],
                   paral=TRUE)

#Build error output
asmmc_error <- 
  mean(c(sum(rf_ms$model.rf$confusion.matrix["1",c("2","3")])/sum(rf_ms$model.rf$confusion.matrix["1",c("1","2","3","6")]),
         sum(rf_ms$model.rf$confusion.matrix["6",c("2","3")])/sum(rf_ms$model.rf$confusion.matrix["6",c("1","2","3","6")])))

modsel1 <- levels(pred_ms$allocation)[as.numeric(pred_ms$allocation)]
levels_rest <- levels(pred_ms$allocation)[-as.numeric(pred_ms$allocation)]
modsel2 <- levels_rest[which.max(pred_ms$vote[-as.numeric(pred_ms$allocation)])]
if (modsel1==6){postmed <- 0;post_lq <- NA;post_hq <- NA} else {
  fittedparam <- switch(modsel1,"1"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
  postmed <- fittedparam$med; post_lq <- fittedparam$quantiles[2]  
  post_hq <- fittedparam$quantiles[5] 
}
if (modsel2==6){postmed2 <- 0} else {
  fitted2param <- switch(modsel2,"1"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
  postmed2 <- fitted2param$med 
}

out1 <- as.data.frame(list(dataset=names(data_main)[i],oob=round(rf_ms$prior.err,3),
                           asmmcoob=round(asmmc_error,3),
                      modelsel=modsel1,postp=round(pred_ms$post.prob,3),postmed = postmed, 
                      postq025=post_lq,
                      postq975 = post_hq,modsel2=modsel2,postmed2=postmed2))

write.table(out1,file=paste0("res/abcres",rep1,"_",names(data_main)[i],".txt"),quote = FALSE,
            row.names = FALSE)




save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/rfres",rep1,names(data_main)[i],".RData"))

