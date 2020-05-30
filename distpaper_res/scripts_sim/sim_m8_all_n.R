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
#source("divfunwrappers.R")

sfsR=function(hapmatrix){
  n=dim(hapmatrix)[1]
  S=dim(hapmatrix)[2]
  I=1:(n-1)
  a1=sum(1/I)
  a2=sum(1/I^2)
  b1=(n+1)/3/(n-1)
  b2=2*(n^2+n+3)/9/n/(n-1)
  c1=b1-1/a1
  c2=b2-(n+2)/a1/n+a2/a1^2
  e1=c1/a1
  e2=c2/(a1^2+a2)
  Dnum=sum(hamming.distance(hapmatrix))/n/(n-1)-S/a1
  Dden=sqrt(e1*S+e2*S*(S-1))
  i=0
  for (j in 1:S){
    i[j]=length(which(hapmatrix[,j]!=hapmatrix[1,j]))
  }
  
  ii=as.numeric(names(table(i)))
  zeta=as.numeric(table(i))
  g1=(n-2)/6/(n-1)
  g2=(18*n^2*(3*n+2)*sum(1/(1:n)^2)-88*n^3-9*n^2+13*n-6)/9/n/(n-1)^2
  theta=S/a1
  theta2=S*(S-1)/(a1^2+b1)
  Hnum=sum(hamming.distance(hapmatrix))/n/(n-1)-sum(ii*zeta)/(n-1)
  Hden=sqrt(g1*theta+g2*theta2)
  
  return(list(TajimaD=Dnum/Dden,FayandWuH=Hnum/Hden))
  
}


require(e1071)

#' Input: 0-1 SNP matrix Output: named vector (TajimaD, FayandWuH)
D_H <- function(seq1){unlist(sfsR(seq1+1))}



divfun_most <- function(seq1,n_ind){
  quantv1 <- seq(.1,.9,.1)
  if (is.matrix(seq1)){
    out1 <- c(quant_hm_oc(seq1,quant_v = quantv1),mean_sd_oc(seq1),
              hammfun(seq1,quant_v = quantv1),phylolength(seq1,quant_v = quantv1),r2fun(seq1,quant_v = quantv1),
              f_nucdiv_S(spectrum01(seq1)),D_H(seq1),allele_freqs(seq1,quant_v = quantv1),spectrum01(seq1))}
  else {out1 <- rep(NA,51+n_ind)}
  names(out1) <-c("hm(O)",paste("O: qu",seq(.1,.9,.1)),"mean(O)","sd(O)",
                  paste("Ham: qu",seq(.1,.9,.1)),
                  paste("Phy: qu",seq(.1,.9,.1)),
                  paste("r2: qu",seq(.1,.9,.1)),"Nucl. div.","S","Tajima's D", "Fay & Wu's H",
                  paste("AF: qu",seq(.1,.9,.1)),paste0("S",1:(n_ind-1)))
  return(out1)}

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
