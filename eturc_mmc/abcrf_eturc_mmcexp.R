#' Test MMC+exp vs best fitting model

#' 2 Args: Which data set from mtb_data_call90.RData, 
#' which replication (= seed set, commented out), how many cores
#' 1st argument "color" code as text string
#' red=BC,green=SC,lightblue=FC,pink=DIV,kenya
#' 2nd argument scaffold 1-18
args1 <- commandArgs(TRUE) 

col1 <- args1[1]
i <- as.integer(args1[2])
log_arg <- TRUE

#' Seeds for replicability - TODO check whether parallelisation causes issues
seed1 <- 1
set.seed(1)
randomseeds <- sample(1:100000,20)
set.seed(randomseeds[1])




library(abcrf)
library(parallel)
library(psych)
library(gap)
library(phyclust)

source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/ext_fun_CREB.R")
source("../general_scripts/misc_scripts.R")
source("../general_scripts/lambdacoal_sim.R")
source("../general_scripts/divstats.R")
source("../distpaper_res/construct_prior.R")
source("../MTB_MMC_repo/serial_sampling_priorsim.R")



#' number of cores for parallel computation
mc1 <- 7#10



gran1 <- c(log(0.5),log(5000))
#gran2 <- c(log(0.5),log(20000))
growthrange0 <- c(log(0.5),log(5000))

betarange0 <- c(1,2)
diracrange0 <- c(0,1)
n0 <- 125000#80000 for the reduced set

mod1 <- 1#df_bestm[i,col1] We always compare w. exp. growth


#beta_gran1 <- c(0,log(5000))
#dirac_gran1 <- c(log(1),log(5000))

beta_gran1 <- c(0,log(50))
dirac_gran1 <- c(log(1),log(50))

#' For Beta(2-a,a), we use the modified Moran scaling (to be consistent when we use full BETA), 
#' see Corollary 1 and Eq 10 Freund 2020. Input is the coalescent parameter a and growth rate g
#' Output is the correct timescale function for Beta coalescents, which is
#' the inverse of (g*a)^(-1)*(exp(g*a*t)-1), since we simulate standard Beta and want the
#' time-changed jump times
#' For Dirac coalescents we just set a=1.5 which corresponds to
#' gamma=1.5 in the Eldon-Wakeley model 
#' see also Matuszewski et al. A different gamma just scales growth
#' differently...

G_tc <- function(g,a){function(t){(a*g)^(-1)*log(a*g*t+1)}}

#' Prior setting for MMC+growth, uses 10K sims to 
#' approximate expected total tree length for Watterson's estimator
#' coal_m="beta" or coal_m="dirac"
#' Model is then numbered 20 for Beta+exp and 30 for Dirac+exp
prior_obs_s_mmc_tc <- function(n_ind,nsimul=100,coal_m="beta",coal_p=c(1,1.9),expg_p=c(1,log(5000)),
                               s_obs = c(10,20),discr=TRUE,log_g=TRUE,mc1=7){
if (length(s_obs)==1){s_obs <- rep(s_obs,2)}
if (coal_m=="beta"){  
theta_f <- function(a,g,s){th1 <- 2*s/est_elength_serial_tc(ct = rep(0,n_ind),
                                                            rate1 = function(n){beta_bc_rates(n,2-a,a)},
                                                            Ginv = G_tc(g,a),mc1 = mc1)}
}
if (coal_m=="dirac"){  
theta_f <- function(p,g,s){th1 <- 2*s/est_elength_serial_tc(ct = rep(0,n_ind),
                                                                rate1 = function(n){dirac_bc_rates(n,p)},
                                                                Ginv = G_tc(g,1.5),mc1 = mc1)}
  }
mod_vec <- rep(switch(coal_m,"beta"=20,"dirac"=30),nsimul)
if (!discr){
sample_coalp <- runif(min = coal_p[1],max = coal_p[2],
                      n = nsimul)
sample_expg <- runif(min = expg_p[1],max = exp_p[2],
                      n = nsimul)
if (log_g){sample_expg <- exp(sample_expg)}
sample_s <- sample(s_obs,nsimul,TRUE)

thef <- function(i){theta_f(sample_coalp[i],sample_expg[i],
                            sample_s[i])}
theta_drawn <- sapply(1:nsimul,thef)

prior1 <- cbind(rep(n_ind,nsimul),mod_vec,sample_coalp,sample_expg,
                sample_s,theta_drawn)
}
if (discr){
temp_g <- seq(expg_p[1],expg_p[2],length.out = 10)
if (log_g){temp_g <- exp(temp_g)}
params_used <- expand.grid(beta=seq(coal_p[1],coal_p[2],length.out = 10),
                           temp_g,
                           s=s_obs)
theta_grid <- sapply(1:nrow(params_used),
                       function(i){theta_f(params_used[i,1],
                                           params_used[i,2],
                                           params_used[i,3])})
params_used <- cbind(params_used,theta_grid)
sample_params <- sample(nrow(params_used),nsimul,TRUE)
prior1 <- cbind(rep(n_ind,nsimul),mod_vec,
                params_used[sample_params,])
}

colnames(prior1) <- c("n_ind","model","coal_param","expg_param",
                        "s_obs","theta_watt")
  return(prior1)
}


#' Wrapper for simulation

sim_seq <- function(nsamp1,theta1,coal_param=0,model=1){
  KM <- (model==1 & coal_param==0) | (model==2 & length(coal_param)==1 & coal_param==2) | (model==6)
  #KM <- KM | (model == 7 & coal_param==0)
  if (KM){raw_seq <- ms(nsam=nsamp1,opts=paste("-t",theta1))
  seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])}
  if (model==1 & !KM & coal_param>=0){
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

#' Diversity statistics

divfun_foldf <- function(seq1){
  if (is.matrix(seq1)){
    out1 <- c(hammfun(seq1),phylolength(seq1),r2fun(seq1),
              f_nucdiv_S(spectrum01(seq1)),
              allele_freqs(seq1,minor=TRUE))} else {
                out1 <- rep(NA,22)}
  names(out1) <-c(paste0("Ham_q",seq(1,9,2)),paste0("Phy_q",seq(1,9,2)),
                  paste0("r2_q",seq(1,9,2)),
                  c("pi","S"),
                  paste0("AF",seq(1,9,2)))
  return(out1)}


load(paste0("../eturc_mmc/input_data/ABCinput_",col1,"_18scaff.RData"))  


#' Get prior
setwd("../general_scripts/")


nsimul_v <- rep(0,6)
nsimul_v[mod1] <- n0 
cat("prior standard \n")
prior1 <- prior_obs_s_cont(n_ind = data_scaffs[[i]]$n_ind,
                      models = mod1,nsimul=nsimul_v,
                      ranges = list(growthrange0,
                                    betarange0,
                                    diracrange0,NULL,NULL,c(0,0)),
                      s_obs = rep(data_scaffs[[i]]$n_mut,2),
                      log_growth = TRUE,
                      include_g0 = 0.02*n0,
                      mc1 = mc1)
prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
cat("prior Betaexp \n")

prior2 <- prior_obs_s_mmc_tc(n_ind = data_scaffs[[i]]$n_ind,
                             nsimul = n0,coal_p = c(1,1.9),coal_m = "beta",
                             expg_p = beta_gran1,
                             s_obs = data_scaffs[[i]]$n_mut,
                             mc1 = mc1,discr = TRUE,log_g = TRUE)
prior2[,"theta_watt"] <- sapply(prior2[,"theta_watt"],log_smear)
rownames(prior2) <- NULL

cat("prior Dirac+exp \n")

prior3 <- prior_obs_s_mmc_tc(n_ind = data_scaffs[[i]]$n_ind,
                             nsimul = n0,coal_p = c(0.05,0.95),
                             coal_m = "dirac",
                             expg_p = dirac_gran1,
                             s_obs = data_scaffs[[i]]$n_mut,
                             mc1 = mc1,discr = TRUE,log_g = TRUE)      
prior3[,"theta_watt"] <- sapply(prior3[,"theta_watt"],log_smear)
rownames(prior3) <- NULL
setwd("../eturc_mmc/")


clu1 <- makeForkCluster(nnodes = mc1)
cat("sim standard \n")

sims1 <- parApply(clu1,prior1,1,function(x){
  divfun_foldf(sim_seq(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],model = x[2]))})

cat("sim Beta+exp \n")

sims1 <- cbind(sims1,parApply(clu1,prior2,1,function(x){
divfun_foldf(Lambda_seq_sim_difft(
           Lambda_difft_coal_tc(coal_times = rep(0,x[1]),
                                ratef1 = function(n){beta_bc_rates(n,2-x[3],x[3])},
                                Ginv = G_tc(g = x[4],a = x[3])),
           theta = x[6]))}))

cat("sim Dirac+exp \n")

sims1 <- cbind(sims1,parApply(clu1,prior3,1,function(x){
    divfun_foldf(Lambda_seq_sim_difft(
    Lambda_difft_coal_tc(coal_times = rep(0,x[1]),
                         ratef1 = function(n){dirac_bc_rates(n,x[3])},
                         Ginv = G_tc(g = x[4],a = 1.5)),
    theta = x[6]))}))
stopCluster(clu1)

prior1 <- rbind(prior1,prior2[,-4])
prior1 <- rbind(prior1,prior3[,-4])

# Save sims for further analysis
save(prior1,prior2,prior3,sims1,file=paste0("sims/sim_mmcg_",col1,i,".RData"))

model1 <- prior1[,2]
coal_p <- prior1[,3]

bad_cols <- apply(sims1,2,function(v){any(is.na(v))})


sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
coal_p <- coal_p[!(bad_cols)]

sumstat1 <- t(sims1)
rownames(sumstat1) <- NULL
target1 <- as.data.frame(t(divfun_foldf(data_scaffs[[i]]$seq_0_1)))
 

model1 <- as.factor(model1)
sumstat1a <- data.frame(model1,t(sims1)) 

good_cols <- lapply(sumstat1a[,-1],
                    function(v){sum(tapply(v,sumstat1a[,1],function(x){range(x)[1]==range(x)[2]}))>1})
good_cols <- !(c(FALSE,unname(unlist(good_cols))))


rf_ms <- abcrf(model1~.,data=sumstat1a[,good_cols],paral = TRUE)

#' Output data  
data2model <- predict(rf_ms,target1[good_cols[-1]],sumstat1a[,good_cols],
                   paral=TRUE)

conf.matrix <- rf_ms$model.rf$confusion.matrix
save(data2model,conf.matrix,
     file=paste0("resdata/abcout_",col1,i,"_mmcg.RData"))

