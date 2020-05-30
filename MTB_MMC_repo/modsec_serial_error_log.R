#' Perform ABC of serial simulations with params mimicking
#' one of three data sets (see below)
#' using a random forest constructed on ultrametric trees 
#' (so non-serial sampling)
#' ABC for model selection
#' Run via Rscript, first argument integer 1-3
#' (1=Eldholm 2015, 2=Lee 2015, 3=Roetzer 2013) 
#' second argument sets prior for growth rates (in [0.5,5000]) 
#' 0=take uniform prior 
#' 1=take uniform prior on the logarithmic scale, i.e.
#'   draw g from range and then use exp(g) as parameter
#'   In particular: Give range as log(range)

args <- commandArgs(TRUE)

#' Params for the data sets mimicked 
load("serial_data.RData")


i <- as.integer(args[1])
log_growth <- as.logical(as.numeric(args[2]))

#' Set seeds (if wanted)
load("seeds_TB.RData")
seed1 <- seeds1[2*23+i]
set.seed(seed1)

#' R packages needed
library(abcrf)
library(parallel)
library(psych)
library(gap)
library(phyclust)

#' Scripts from the repo
source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/ext_fun_CREB.R")
source("../general_scripts/misc_scripts.R")
source("../general_scripts/lambdacoal_sim.R")
source("../general_scripts/divstats.R")
source("../distpaper_res/construct_prior.R")

#' How many cores used for parallel computation?
mc1 <- 7#16

#' Parameter range for prior distributions for simulating
#' ultrametric coalescent trees from three model classes
#' Kingman+exp, Beta, Dirac, from which the random forest is built
if (log_growth){growthrange0 <- c(log(0.5),log(5000))} else {#seq(0,2500,1)
growthrange0 <- c(0.5,5000)}
betarange0 <- c(1,2)#seq(1,1.975,0.025)
diracrange0 <- c(0,1)#seq(0.025,0.975,0.025)
n0 <- 125000 #Number of coalescent simulations per model class


#' Build prior distribution in all model classes (including standard Kingman)
#' model 1: Kingman+exp, 2: Beta, 3: Dirac, 6: Kingman
setwd("../general_scripts/")
prior_rf <- prior_obs_s_cont(params_serial[[i]]$n_ind,
                             models=c(1,2,3,6),
                             nsimul=c(n0,n0,n0,0,0,n0),
                            ranges = list(growthrange0,
                                    betarange0,diracrange0,NULL,NULL,c(0,0)),
                      s_obs = rep(params_serial[[i]]$n_mut,2),
                      log_growth = log_growth,
                      include_g0 = 0,mc1 = mc1)

prior_rf[,"theta_watt"] <- sapply(prior_rf[,"theta_watt"],log_smear)
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

#' Parallel simulation (non-serial)
clu1 <- makeForkCluster(nnodes = mc1)
sims_rf <- parApply(clu1,prior_rf,1,function(x){
  divfun_ohaf(sim_seq(nsamp1 = x[1],theta1 = x[5],
                       coal_param = x[3],model = x[2]))})
stopCluster(clu1)


model_rf <- prior_rf[,2]
coal_rf <- prior_rf[,3]

#' For safety: Remove simulations with NA (due to overflow)
bad_cols <- apply(sims_rf,2,function(v){any(is.na(v))})
sims_rf <- sims_rf[,!(bad_cols)]
model_rf <- model_rf[!(bad_cols)]



model_rf <- as.factor(model_rf)
sumstat_rf <- data.frame(model_rf,t(sims_rf)) 


good_cols <- lapply(sumstat_rf[,-1],
                    function(v){sum(tapply(v,sumstat_rf[,1],function(x){range(x)[1]==range(x)[2]}))>1})
good_cols <- !(c(FALSE,unname(unlist(good_cols))))

#' Construct random forest (same random forest for all serial sim parameters)
rf_ms <- abcrf(model_rf~.,data=sumstat_rf[,good_cols],paral = TRUE,ncores = mc1)

#' Load pre-simulated simulations under serial coalescent models
load(switch(i,"resdata/serial_simsEldholm2015.RData",
            "resdata/serial_simsLee2015.RData",
            "resdata/serial_simsRoetzer2013.RData"))



modsel_rf_serial <- vector("list",length(sims1))
names(modsel_rf_serial) <- names(sims1)

err_rf_serial <- modsel_rf_serial
errmmc_rf_serial <- modsel_rf_serial

# For each parameter combination of serial coalescent simulation,
# perform ABC model selection 
for (k in 1:length(sims1)){
#' Record misclassification both as any other model (not K+exp as Kingman, 
#' since model classes overlap)
checkmod <- switch(k,c("2","3","6"),c("1","3","6"),c("1","2","6"))
checkmmc <- switch(k,c("2","3"),c("1","6"),c("1","6"))
  
modsel_rf_serial[[k]] <- vector("list",nrow(params_order[[k]]))  
err_rf_serial[[k]] <- matrix(-1,ncol=length(unique(params_order[[k]][,1])),
                                nrow=length(unique(params_order[[k]][,2])))
colnames(err_rf_serial[[k]]) <- sort(unique(params_order[[k]][,1])) #cols: coalescent param
rownames(err_rf_serial[[k]]) <- sort(unique(params_order[[k]][,2])) #rows: scaling factor
errmmc_rf_serial[[k]] <- err_rf_serial[[k]] 

for (j in 1:nrow(params_order[[k]])){
#' For safety: Remove simulations with NA (due to overflow)
  bad_cols <- apply(sims1[[k]][[j]],2,function(v){any(is.na(v))})
bad_cols <- bad_cols 
  
prior_mtest <- prior1[[k]][[j]][!bad_cols,]
sims_mtest <- sims1[[k]][[j]][,!(bad_cols)]
model_mtest <- as.factor(prior_mtest[,2])
  
stats_mtest <- data.frame(model_mtest,t(sims_mtest))
  
#' ABC via random forest: model selection
pred_m_rf <- predict(rf_ms,stats_mtest[,good_cols],sumstat_rf[,good_cols],paral = TRUE,
                     ncores = mc1)

#' Record errors
modsel_rf_serial[[k]][[j]] <- pred_m_rf

err_rf_serial[[k]][as.character(params_order[[k]][j,2]),
                   as.character(params_order[[k]][j,1])] <-  
sum(pred_m_rf$allocation %in% checkmod)/length(pred_m_rf$allocation)

errmmc_rf_serial[[k]][as.character(params_order[[k]][j,2]),
                   as.character(params_order[[k]][j,1])] <-  
sum(pred_m_rf$allocation %in% checkmmc)/length(pred_m_rf$allocation)
}  
}

if (log_growth){addon <- "logg"} else {addon <- "contg"}

save(errmmc_rf_serial,err_rf_serial,modsel_rf_serial,
     file=paste0("resdata/serial_modsec_",names(params_serial)[i],addon,".RData"))

