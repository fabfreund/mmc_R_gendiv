#' Simulate 1000 samples matching to three serial sampled data sets 
#' under serial sampling with sampling times within some  
#' percentage of total length of ultrametric tree and fixed growth rate 


args <- commandArgs(TRUE) # Data set nr., but from the serial ones (so 1:3)
datasetnr <- as.numeric(args[1]) 

#' R packages needed
library(parallel)
library(psych)
library(extraDistr)

#' R code needed
source("../general_scripts/lambdacoal_sim.R")
source("../general_scripts/divstats.R")
source("serial_sampling_priorsim.R") #Main tool, this is the serial simulation
                                     #script

#' Load parameters from the three data sets Lee 2015, Eldholm 2015, Roetzer 2013
#' (we need sample size, number of mutations, real sampling times)
load("serial_data.RData")

#' Diversity function used 



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

#' Params for time scaling and coalescent models
coal_param <- list(c(1,10,50,100,250,500,1000,2000), #growth rates
                   seq(1,2,.2), #alpha for Beta coalescent, includes Kingman
                   seq(.1,.9,.2)) #Dirac coalescent param



c1 <- c(seq(0,0.5,.1),.75,1,1.5)  #coalescent time proportion that is used for the 
                                  #serial sampling
kexp1 <- c(TRUE,FALSE,FALSE) #Only growth for kexp

#' Output files
sims1 <- list(vector("list",length(coal_param)*length(c1)),
              vector("list",length(coal_param)*length(c1)),
              vector("list",length(coal_param)*length(c1)))
prior1 <- list(vector("list",length(coal_param)*length(c1)),
               vector("list",length(coal_param)*length(c1)),
               vector("list",length(coal_param)*length(c1)))
names(sims1) <- c("K+exp","Beta","Dirac")
names(prior1) <- c("K+exp","Beta","Dirac")

#' Simulation parameters
n0 <- 1000 #Check 1K sims per condition cell whether it gets correctly sorted
mc1 <- 14#16


#' Extract data parameters
n_mut <- params_serial[[datasetnr]]$n_mut
n1 <- params_serial[[datasetnr]]$n_ind
sampl_t <- unname(params_serial[[datasetnr]]$sample_t) #Time back from sampling date,coalescent wants them positive 
sampl_t <- sampl_t/max(sampl_t)



for (mod1 in 1:3){ #Run through 3 coalescent models

k <- 0
for (g2 in coal_param[[mod1]]){
  for (c2 in c1){
k <- k+1
cat("Params",g2,c2,"\n")
coal_rate_f <- switch(mod1,function(n){beta_bc_rates(n,0,2)},
                      function(n){beta_bc_rates(n,2-g2,g2)},
                      function(n){dirac_bc_rates(n,g2)})
coal_h <- est_eheight_serial(ct = rep(0,length(sampl_t)),
                             rate1 = coal_rate_f,kme = kexp1[mod1],
                             rho = g2,mc1 = mc1)
coal_l <- est_elength_serial(c2*coal_h*sampl_t,coal_rate_f,kexp1[mod1],
                             g2,mc1 = mc1)  
watt_est <- 2*n_mut/coal_l


# x[1] is g2, x[2] is c2
coal_t_exp <- function(x){
  model_bool <- TRUE
  temp1 <- Lambda_difft_coal(coal_times = x[2]*coal_h*sampl_t,
                             ratef1 = coal_rate_f,
                             kmexp = kexp1[mod1],rho = x[1])
  return(Lambda_seq_sim_difft(temp1,watt_est))
} 

mod_numb <- switch(mod1,9,10,11) 
#Model 9: serial growth, 10 serial beta, 11 serial Dirac

prior1[[mod1]][[k]] <- cbind("model"=rep(mod_numb,n0),
                             "coal_param"=rep(g2,n0),
                             "coal_scale"=rep(c2,n0))
clu1 <- makeForkCluster(nnodes = mc1)
sims1[[mod1]][[k]] <- parApply(clu1,prior1[[mod1]][[k]][,-1],1,
                               function(x){divfun_ohaf(coal_t_exp(x))})
stopCluster(clu1)
  }}
}

params_order <- vector("list",3)
names(params_order) <- c("K+exp","Beta","Dirac")
for (mod1 in 1:3){
params_order[[mod1]] <- expand.grid(c1,coal_param[[mod1]])
params_order[[mod1]] <- params_order[[mod1]][,c(2,1)]
colnames(params_order[[mod1]]) <- c("coal_param","coal_scale")
}

save(prior1,sims1,params_order,
     file=paste0("resdata/serial_sims",names(params_serial)[datasetnr],".RData"))

