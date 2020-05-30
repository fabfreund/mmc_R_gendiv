#' Performs ABC model selection KM+exp w. pop. structure vs. BETA vs. Dirac
#' for Lee 2015

load("mtb_data_call90.RData")

#' We only need it for data set Lee 2015

i <- which(names(data_main)=="Lee2015")


#' Seeds

rep1 <- 2

#' Seed setting, only partially applied in the manuscript
load("seeds_TB.RData")
seed1 <- seeds1[29*4+i+rep1] 
set.seed(seed1)

#' We perform a random forest based ABC for model selection and parameter estimation
#' based on minimal observable clade sizes, Hamming distances and AF+
#' We perform model selection between Beta, Dirac and structured Kingman + exp

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

#' number of cores for parallel computation
mc1 <- 7


gran1 <- seq(0,5000,2)
gran2 <- seq(0,20000,5)
f1 <- function(x){switch(x,gran1,gran2)}
growthrange0 <- lapply(c(1,2,1,1,2,1,2,2,1,1,1,2,2,1,1,1,2,2,2,1,1,1,1,1,2,2,2,2,2),f1)
betarange0 <- c(1,2)
diracrange0 <- c(0,1)
migrange <- c(.25,.5,1,2,3)
n0 <- 125000
n_1mig <- n0/length(migrange) 


#' Prior for K+exp+struct, mutation set as Watterson's estimate
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
  coal_p <- ranges[sample_g]
  prior1 <- cbind(rep(n_ind,sum(nsimul)),mod_vec,coal_p,
                  s_obs[sample_s],sample_theta)
  colnames(prior1) <- c("n_ind","model","coal_param","s_obs","theta_watt")
  return(prior1)
}



setwd("../general_scripts/")
prior1 <- prior_obs_s_cont(data_main[[i]]$n_ind,models=c(2,3),nsimul=c(0,n0,n0,0,0,0),
                      ranges = list(growthrange0[[i]],
                                    betarange0,diracrange0,NULL,NULL,c(0,0)),
                      s_obs = rep(data_main[[i]]$n_mut,2),mc1 = mc1)

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)

prior_temp <- vector("list",length(migrange))
prior_temp2 <- NULL
l <- 0
for (m1 in migrange){
  l <- l+1
  prior_temp[[l]]  <- prior_obs_s_exp_popstr(n_ind = data_main[[i]]$n_ind,nsimul = n_1mig,
                                             s_obs = rep(data_main[[i]]$n_mut,2),migr = m1,samples = c(61,36,49,1))
  
  prior_temp[[l]][,"theta_watt"] <-  sapply(prior_temp[[l]][,"theta_watt"],log_smear)
  prior_temp2 <- rbind(prior_temp2,prior_temp[[l]])
}


setwd("../MTB_MMC_repo/")

#' Simulation for Beta & Dirac

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

#' Simulate K+exp+struct

sim_seq_popstr <- function(nsamp1=100,theta1,coal_param=0.5,samples=c(50,50),migr=0.5){
  ms_com <- paste("-t",theta1,"-G",coal_param,"-I",length(samples),paste(samples,collapse = " "),migr)
  raw_seq <- ms(nsam=nsamp1,opts=ms_com)
  seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])
  if (sum(seq1)==0){return(NA)} 
  return(seq1)  
}

#' Diversity stats

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

sims2 <- NULL

clu1 <- makeForkCluster(nnodes = mc1)

sims1 <- parApply(clu1,prior1,1,function(x){
  divfun_ohaf(sim_seq(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],model = x[2]))})

for (l in 1:length(migrange)){
  sims2 <- cbind(sims2,parApply(clu1,prior_temp[[l]],1,function(x){
    divfun_ohaf(sim_seq_popstr(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],samples = c(61,36,49,1),
                                 migr = migrange[l]))}))  
}


stopCluster(clu1)

prior1 <- rbind(prior1,prior_temp2)

sims1 <- cbind(sims1,sims2)


model1 <- prior1[,2]
coal_p <- prior1[,3]

bad_cols <- apply(sims1,2,function(v){any(is.na(v))})


sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
coal_p <- coal_p[!(bad_cols)]

sumstat1 <- t(sims1)
rownames(sumstat1) <- NULL
target1 <- as.data.frame(t(divfun_ohaf(data_main[[i]]$seq_0_1)))
 
#quantiles to be returned for posterior distribution
quant_ret <- seq(0,1,0.025)
#' K+exp+struct is model 8, we just fit growth, mig is nuisance param
param1 <- coal_p[model1==8]
sumstat2 <- data.frame(param1,sumstat1[model1==8,])
growth_rf <- regAbcrf(param1~.,
                      data=sumstat2,
                      paral = TRUE)
growthfit_rf <- predict(growth_rf,target1,
                        sumstat2,paral=TRUE,
                        quantiles = quant_ret)

print(growthfit_rf)

param1 <- coal_p[model1==2]
sumstat2 <- data.frame(param1,sumstat1[model1==2,])  
beta_rf <- regAbcrf(param1~.,
                    data=sumstat2,
                    paral = TRUE)
betafit_rf <- predict(beta_rf,target1,sumstat2,paral=TRUE,
                      quantiles = quant_ret)

print(betafit_rf)

param1 <- coal_p[model1==3]
sumstat2 <- data.frame(param1,sumstat1[model1==3,])   
dirac_rf <- regAbcrf(param1~.,
                     data=sumstat2,
                     paral = TRUE)
diracfit_rf <- predict(dirac_rf,target1,sumstat2,
                       paral=TRUE,quantiles = quant_ret)

print(diracfit_rf)

model1 <- as.factor(model1)
sumstat1a <- data.frame(model1,t(sims1)) 

good_cols <- lapply(sumstat1a[,-1],
                    function(v){sum(tapply(v,sumstat1a[,1],function(x){range(x)[1]==range(x)[2]}))>1})
good_cols <- !(c(FALSE,unname(unlist(good_cols))))


rf_ms <- abcrf(model1~.,data=sumstat1a[,good_cols],paral = TRUE)


print(rf_ms)
  
pred_ms <- predict(rf_ms,target1[good_cols[-1]],sumstat1a[,good_cols],
                   paral=TRUE)

print(pred_ms)


asmmc_error <- sum(rf_ms$model.rf$confusion.matrix["8",
                    c("2","3")])/sum(rf_ms$model.rf$confusion.matrix["8",c("2","3","8")])


modsel1 <- levels(pred_ms$allocation)[as.numeric(pred_ms$allocation)]
levels_rest <- levels(pred_ms$allocation)[-as.numeric(pred_ms$allocation)]
modsel2 <- levels_rest[which.max(pred_ms$vote[-as.numeric(pred_ms$allocation)])]
fittedparam <- switch(modsel1,"8"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
postmed <- fittedparam$med; post_q <- fittedparam$quantiles

  fitted2param <- switch(modsel2,"8"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
  postmed2 <- fitted2param$med 

out1 <- as.data.frame(as.list(c(dataset=names(data_main)[i],oob=round(rf_ms$prior.err,3),
                           asmmcoob=round(asmmc_error,3),
                      modelsel=modsel1,postp=round(pred_ms$post.prob,3),postmed = postmed, 
                      qu=post_q,modsel2=modsel2,postmed2=postmed2)))


write.table(out1,file=paste0("res/abcres",rep1,"_",names(data_main)[i],"_struct.txt"),quote = FALSE,
            row.names = FALSE)


