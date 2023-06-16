#' ABC-RF on samples from different Exserohilum turcicum populations
#' 1st argument "color" code as text string
#' red=BC,green=SC,lightblue=FC,pink=DIV,kenya
#' 2nd argument scaffold 1-18
#' 3rd argument T,F: log prior TRUE or FALSE
args1 <- commandArgs(TRUE) 

col1 <- args1[1]
i <- as.integer(args1[2])
log_arg <- as.logical(args1[3])

setwd("../MTB_MMC_repo/")



seed1 <- 1
set.seed(1)
randomseeds <- sample(1:100000,20)
set.seed(randomseeds[7])
#set.seed(randomseeds[5])
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


#' number of cores for parallel computation
mc1 <- 3#10


#' Simulation parameter for model classes for the ABC
if (log_arg){gran1 <- c(log(0.5),log(5000))} else {gran1 <- c(0,5000)}
growthrange0 <- gran1
betarange0 <- c(1,2)
diracrange0 <- c(0,1)
n0 <- 125000


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
divfun_foldf <- function(seq1){
  if (is.matrix(seq1)){
    out1 <- c(hammfun(seq1),phylolength(seq1),r2fun(seq1),
              f_nucdiv_S(spectrum01(seq1)),
              allele_freqs(seq1,minor=TRUE))} else {
                out1 <- rep(NA,22)}
  names(out1) <-c(paste0("Ham_q",seq(1,9,2)),paste0("Phy_q",seq(1,9,2)),
                  paste0("r2_q",seq(1,9,2)),
                  c("pi","S"),
                  paste0("MAF_q",seq(1,9,2)))
  return(out1)}


#' Load data
load(paste0("../eturc_mmc/input_data/ABCinput_",col1,
            "_18scaff.RData"))  

#' Simulate
setwd("../general_scripts/")

prior1 <- prior_obs_s_cont2(n_ind = data_scaffs[[i]]$n_ind,
                            models=c(1,2,3),
                            nsimul=c(n0,n0,n0,0,0,0),
                            ranges = list(growthrange0,
                                          betarange0,diracrange0,NULL,NULL,c(0,0)),
                            s_obs = rep(data_scaffs[[i]]$n_mut,2),
                            log_growth = log_arg,include_g0 = 0.02*n0,
                            include_bsz=n0*0.05,
                            mc1 = mc1)

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
setwd("../eturc_mmc/")


#' Parallelized coalescent simulations
clu1 <- makeForkCluster(nnodes = mc1)
sims1 <- parApply(clu1,prior1,1,function(x){
  divfun_foldf(sim_seq(nsamp1 = x[1],theta1 = x[5],
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
target1 <- as.data.frame(t(divfun_foldf(data_scaffs[[i]]$seq_0_1)))

#quantiles to be returned for posterior distribution
quant_ret <- seq(0,1,0.025) 

#cleanup
#save.image("test.RData") #for internal testing
gc()


#' Parameter estimation Kingman+exp
param1 <- coal_p[model1==1]
param_g <- param1
sumstat2 <- data.frame(param1,sumstat1[model1==1,])
growth_rf <- regAbcrf(param1~.,
                      data=sumstat2,
                      paral = TRUE,ncores=mc1)
growthfit_rf <- predict(growth_rf,target1,
                        sumstat2,paral=TRUE,ncores=mc1,
                        quantiles = quant_ret)

#cleanup
gc()

#' Parameter estimation Beta
param1 <- coal_p[model1==2]
sumstat2 <- data.frame(param1,sumstat1[model1==2,])  
beta_rf <- regAbcrf(param1~.,
                    data=sumstat2,
                    paral = TRUE,ncores=mc1)
betafit_rf <- predict(beta_rf,target1,sumstat2,paral=TRUE,ncores=mc1,
                      quantiles = quant_ret)

#cleanup
gc()

#' Parameter estimation Dirac
param1 <- coal_p[model1==3]
sumstat2 <- data.frame(param1,sumstat1[model1==3,])   
dirac_rf <- regAbcrf(param1~.,
                     data=sumstat2,
                     paral = TRUE,ncores=mc1)
diracfit_rf <- predict(dirac_rf,target1,sumstat2,
                       paral=TRUE,ncores=mc1,quantiles = quant_ret)

#Record errors

#' Compute nmae without 0 values...
compute_no0_nmae <- function(rf,param1){
  which0 <- (param1==0)
  estim1 <- rf$model.rf$predictions
  #next line: should /r be within brackets? 
  nmae_no0 <- mean(abs((estim1[!which0]-param1[!which0])/param1[!which0]))
  mae0 <- mean(abs(estim1[which0]-param1[which0]))
  return(c("nmae_no0"=nmae_no0,"mae0"=mae0))
}

#Make error tables
out2 <- as.data.frame(as.list(c(
  compute_no0_nmae(rf = growth_rf,param1 = param_g),   
  growth_predicterror=growth_rf$model.rf$prediction.error,
  growth_nmae_old=growth_rf$model.rf$NMAE,
  growth_r2=growth_rf$model.rf$r.squared,
  beta_predicterror=beta_rf$model.rf$prediction.error,
  beta_nmae=beta_rf$model.rf$NMAE,
  beta_r2=beta_rf$model.rf$r.squared,
  dirac_predicterror=dirac_rf$model.rf$prediction.error,
  dirac_nmae=dirac_rf$model.rf$NMAE,
  dirac_r2=dirac_rf$model.rf$r.squared)))

write.table(out2,file=paste0("../eturc_mmc/priorerr/err_",names(data_scaffs)[i],
                             ".txt"),quote = FALSE,row.names = FALSE)


#cleanup
rm(dirac_rf,beta_rf,growth_rf)
gc()

model1 <- as.factor(model1)
sumstat1a <- data.frame(model1,t(sims1)) 
#' Remove summary statistics with no variation within a model class
good_cols <- lapply(sumstat1a[,-1],
                    function(v){sum(tapply(v,sumstat1a[,1],function(x){range(x)[1]==range(x)[2]}))>=1}) #added that single no var classes already cause problems...
good_cols <- !(c(FALSE,unname(unlist(good_cols))))

#' ABC model selection
rf_ms <- abcrf(model1~.,data=sumstat1a[,good_cols],paral = TRUE,ncores=mc1)
gc()
pred_ms <- predict(rf_ms,target1[(good_cols[-1])],sumstat1a[,good_cols],
                   paral=TRUE,ncores=mc1)
gc()
#Build error output
asmmc_error <- sum(rf_ms$model.rf$confusion.matrix["1",c("2","3")])/sum(rf_ms$model.rf$confusion.matrix["1",c("1","2","3")])

askm_error <- 
  mean(c(sum(rf_ms$model.rf$confusion.matrix["2","1"])/sum(rf_ms$model.rf$confusion.matrix["2",c("1","2","3")]),
         sum(rf_ms$model.rf$confusion.matrix["3","1"])/sum(rf_ms$model.rf$confusion.matrix["3",c("1","2","3")])))


modsel1 <- levels(pred_ms$allocation)[as.numeric(pred_ms$allocation)]
levels_rest <- levels(pred_ms$allocation)[-as.numeric(pred_ms$allocation)]
modsel2 <- levels_rest[which.max(pred_ms$vote[-as.numeric(pred_ms$allocation)])]
if (modsel1==6){postmed <- 0;post_q <- rep(NA,length(quant_ret))} else {
  fittedparam <- switch(modsel1,"1"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
  postmed <- fittedparam$med; post_q <- fittedparam$quantiles  
}
if (modsel2==6){postmed2 <- 0} else {
  fitted2param <- switch(modsel2,"1"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
  postmed2 <- fitted2param$med 
}

out1 <- as.data.frame(as.list(c(dataset=names(data_scaffs)[i],oob=round(rf_ms$prior.err,3),
                           asmmcoob=round(asmmc_error,3),askmoob=round(askm_error,3),
                      modelsel=modsel1,postp=round(pred_ms$post.prob,3),
                      postmed = postmed, 
                      qu=post_q,modsel2=modsel2,postmed2=postmed2)))

if (log_arg){
write.table(out1,file=paste0("../eturc_mmc/res/res_",names(data_scaffs)[i],
                             ".txt"),quote = FALSE,row.names = FALSE)} else {
write.table(out1,file=paste0("../eturc_mmc/res2/res_",names(data_scaffs)[i],
                             ".txt"),quote = FALSE,row.names = FALSE)}

