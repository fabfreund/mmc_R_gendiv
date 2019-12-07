setwd("../") #To follow the folder structure in the repo

#' We run via Rscript with an argument file which contains the 
#' population samples, migration rates and the seed
#' This then gets loaded
args1 <- commandArgs(TRUE)
source(paste0("scripts_sim/",args1[1])) #To follow the folder structure



#' Source scripts

source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/betaxicoal_sim.R")
source("../general_scripts/divstats.R")
source("construct_prior.R")
source("divfunwrappers.R")

library(parallel)
library(phyclust)
library(gap)

#' Prior function for K+exp+struct, mutation set as Watterson's estimate
#' We approximate the expected total length via simulation of 10K trees
#' 
prior_obs_s_exp_popstr <- function(n_ind,nsimul=100,
                                   ranges = c(0.5,1,2.5,4,7,10,25,50,75,100,500,1000),
                                   s_obs = c(75,75),migr = 0.5,samples=c(50,50)){
  watt_popstruct <- function(g){
    ms_com <- paste("-t 0.000001 -L -G",g,"-I",length(samples),paste(samples,collapse = " "),
                    migr)
    seq1 <- read.ms.output(ms(n_ind,10000,ms_com),is.file = FALSE)
    return(2*mean(seq1$times[,2]))} #ms uses 4N scales instead of 2N scale
  theta <- sapply(ranges,function(g){th1 <- 2*s_obs/watt_popstruct(g)
  return(th1)})  
  sample_g <- sample(seq(along=ranges),nsimul,replace = TRUE)
  sample_s <- sample(seq(along=s_obs),nsimul,replace = TRUE)
  
  sample_theta <- rep(0,nsimul)
  for (i in 1:nsimul){sample_theta[i] <- theta[sample_s[i],sample_g[i]]}
  
  mod_vec <- rep(8,nsimul)
  mig_vec <- rep(migr,nsimul)
  coal_p <- ranges[sample_g]
  popstruct_matrix <- matrix(samples,nrow=nsimul,ncol=length(samples),byrow = TRUE)
  prior1 <- cbind(rep(n_ind,sum(nsimul)),mod_vec,coal_p,
                  s_obs[sample_s],sample_theta,mig_vec,popstruct_matrix)
  colnames(prior1) <- c("n_ind","model","coal_param","s_obs","theta_watt",
                        "migrate",paste("subpop",1:length(samples)))
  return(prior1)
}

#' Simulation function K+exp+struct

sim_seq_popstr <- function(nsamp1=100,theta1,coal_param=0.5,samples=c(50,50),migr=0.5){
  ms_com <- paste("-t",theta1,"-G",coal_param,"-I",length(samples),paste(samples,collapse = " "),migr)
  raw_seq <- ms(nsam=nsamp1,opts=ms_com)
  seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])
  if (sum(seq1)==0){return(NA)} 
  return(seq1)  
}

#' Diversity stats divfun_most


for (i in 1:2){
prior_full <- NULL
for (mig2 in mig1){
   for (popstruct2 in popstruct){
prior_full <- rbind(prior_full,prior_obs_s_exp_popstr(n_ind = nsamp,nsimul = nsim,
                                              ranges = expg_range,
                                              s_obs = s_obs1,migr = mig2,
                                              samples = popstruct2))
                                 }}

clu1 <- makeForkCluster(nnodes = mc1)

sims1 <- parApply(clu1,prior_full,1,function(x){
  divfun_most(sim_seq_popstr(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],
                             samples = x[7:length(x)],migr = x[6]),n_ind = x[1])})

stopCluster(clu1)

prior1 <- prior_full[,1:5]
save(prior1,prior_full,sims1,file=paste0("sims_rep",i,"/",simname))
}
