#' 1st arg data set number, 2nd arg seed set,
#' 3rd arg number of cores for parallelization 
#' 4th arg whether ABC for main analysis should 
#' be directly conducted
args <- commandArgs(TRUE)
i <- as.integer(args[1])
rep1 <- as.integer(args[2])
mc1 <- as.integer(args[3])


load("seeds_TB.RData")
seed1 <- seeds1[23*(rep1-1)+i]
set.seed(seed1)

#' R packages 
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
#' statistics to compute
divfun_ohaf <- function(seq1,private=TRUE){
  if (is.matrix(seq1)){
    log1 <- rep(TRUE,ncol(seq1))
    if (!private){log1 <- (colSums(seq1)>1)}
    if (sum(log1)==0){out1 <- rep(NA,24)} else {
      seq1 <- as.matrix(seq1[,log1])
      out1 <- c(quant_hm_oc(seq1),mean_sd_oc(seq1),hammfun(seq1),
                f_nucdiv_S(spectrum01(seq1)),allele_freqs(seq1,
                                                          quant_v = seq(.1,.9,.1)))}}
  else {out1 <- rep(NA,24)}
  names(out1) <-c("o_hm",paste0("o_q",seq(1,9,2)),"o_mean","o_sd",
                  paste0("ham_q",seq(1,9,2)),
                  "nucdiv","S",
                  paste0("AF_q",seq(1,9,1)))
  return(out1)}
#' simulation script 
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

choose_gr <- function(dataset1,log_growth){
if (!log_growth){
  gran1 <- c(0.5,5000)
  gran2 <- c(0.5,20000)} else {
    gran1 <- c(log(0.5),log(5000))
    gran2 <- c(log(0.5),log(20000))
  } 
f1 <- function(x){switch(x,gran1,gran2)}
growthrange0 <- lapply(c(1,2,2,1,1,1,1,2,1,1,1,1,2,1,2,1,1,1,1,1,
                         1,2,1),f1)
return(growthrange0[[dataset1]])}



n0 <- 125000 # Number of simulations per model class

#' Load the 90% TB data
load("TB_data_ABC_data_fasta.RData")

prior_90 <- vector("list",7)
sims_90 <- vector("list",7)
names(prior_90) <- c("m1_unif","m2_12","m3","m6",
                     "m1_log","m2_02",
                     "m2_02bsz")
names(sims_90) <- names(prior_90)
#' KM+exp w. uniform prior
setwd("../general_scripts/")
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(1),
                           nsimul=c(n0,0,0,0,0,0),
                           ranges = list(choose_gr(i,FALSE),
                                         NULL,NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = FALSE,include_g0 = 0,
                           mc1 = mc1)

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m1_unif"]] <- prior1

#' Beta w. alpha in cont. [1,2] 
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(2),
                           nsimul=c(0,n0,0,0,0,0),
                           ranges = list(NULL,
                                         c(1,2),NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m2_12"]] <- prior1

#' Dirac 
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(3),
                           nsimul=c(0,0,n0,0,0,0),
                           ranges = list(NULL,
                                         NULL,c(0,1),NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m3"]] <- prior1
#' KM
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(6),
                           nsimul=c(0,0,0,0,0,n0),
                           ranges = list(NULL,
                                         NULL,NULL,NULL,NULL,c(0,0)),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m6"]] <- prior1
setwd("../MTB_MMC_repo/")

#' Simulation (should be seeded correctly)
for (j in 1:4){
  name1 <-names(prior_90)[j] 
  clu1 <- makeForkCluster(nnodes = mc1)
  sims_90[[name1]] <- parApply(clu1,prior_90[[name1]],1,function(x){
    divfun_ohaf(sim_seq(nsamp1 = x[1],theta1 = x[5],
                        coal_param = x[3],model = x[2]))})
  stopCluster(clu1)}

sim_seed <- .Random.seed

#' Continue priors
setwd("../general_scripts/")
#' KM+exp w. log prior
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(1),
                           nsimul=c(n0,0,0,0,0,0),
                           ranges = list(choose_gr(i,TRUE),
                                         NULL,NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m1_log"]] <- prior1

#' Beta w. alpha in cont. [0,2] 
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(2),
                           nsimul=c(0,n0,0,0,0,0),
                           ranges = list(NULL,
                                         c(0,2),NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m2_02"]] <- prior1

#' Beta w. alpha cont in [0,2] + 1% BSZ atom
prior1 <- prior_obs_s_cont2(n_ind = data_main[[i]]$n_ind,
                            models=c(2),
                            nsimul=c(0,n0,0,0,0,0),
                            ranges = list(NULL,
                                          c(0,2),NULL,NULL,NULL,NULL),
                            s_obs = rep(data_main[[i]]$n_mut,2),
                            log_growth = TRUE,include_g0 = 0,
                            include_bsz=n0*0.01,
                            mc1 = mc1)

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_90[["m2_02bsz"]] <- prior1

setwd("../MTB_MMC_repo/")

#' simulation

for (j in 5:7){
name1 <-names(prior_90)[j] 
clu1 <- makeForkCluster(nnodes = mc1)
sims_90[[name1]] <- parApply(clu1,prior_90[[name1]],1,function(x){
  divfun_ohaf(sim_seq(nsamp1 = x[1],theta1 = x[5],
                 coal_param = x[3],model = x[2]))})
stopCluster(clu1)}

#' Load the 75% TB data
load("TB_data_ABC_data_het75.RData")
prior_75 <- vector("list",6)
sims_75 <- vector("list",6)
names(prior_75) <- c("m1_unif","m2_12",
                     "m3","m6","m1_log","m2_02")
names(sims_75) <- names(prior_75)
#' KM+exp w. uniform prior
setwd("../general_scripts/")
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(1),
                           nsimul=c(n0,0,0,0,0,0),
                           ranges = list(choose_gr(i,FALSE),
                                         NULL,NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = FALSE,include_g0 = 0,
                           mc1 = mc1)

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_75[["m1_unif"]] <- prior1
#' KM+exp w. log prior
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(1),
                           nsimul=c(n0,0,0,0,0,0),
                           ranges = list(choose_gr(i,TRUE),
                                         NULL,NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_75[["m1_log"]] <- prior1
#' Beta w. alpha in cont. [1,2] 
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(2),
                           nsimul=c(0,n0,0,0,0,0),
                           ranges = list(NULL,
                                         c(1,2),NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_75[["m2_12"]] <- prior1
#' Beta w. alpha in cont. [0,2] 
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(2),
                           nsimul=c(0,n0,0,0,0,0),
                           ranges = list(NULL,
                                         c(0,2),NULL,NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_75[["m2_02"]] <- prior1
#' Dirac 
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(3),
                           nsimul=c(0,0,n0,0,0,0),
                           ranges = list(NULL,
                                         NULL,c(0,1),NULL,NULL,NULL),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_75[["m3"]] <- prior1
#' KM
prior1 <- prior_obs_s_cont(n_ind = data_main[[i]]$n_ind,
                           models=c(6),
                           nsimul=c(0,0,0,0,0,n0),
                           ranges = list(NULL,
                                         NULL,NULL,NULL,NULL,c(0,0)),
                           s_obs = rep(data_main[[i]]$n_mut,2),
                           log_growth = TRUE,include_g0 = 0,
                           mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
prior_75[["m6"]] <- prior1
setwd("../MTB_MMC_repo/")

for (j in seq(along=prior_75)){
  name1 <-names(prior_75)[j] 
clu1 <- makeForkCluster(nnodes = mc1)
sims_75[[name1]] <- parApply(clu1,prior_75[[name1]],1,function(x){
  divfun_ohaf(sim_seq(nsamp1 = x[1],theta1 = x[5],
                      coal_param = x[3],model = x[2]))})
stopCluster(clu1)
}
save(prior_75,prior_90,sims_75,sims_90,sim_seed,
     file=paste0("sims/fullsims_rep",rep1,"_dataset",i,".RData"))