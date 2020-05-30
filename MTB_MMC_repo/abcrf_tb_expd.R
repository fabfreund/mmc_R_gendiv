#' 2nd ABC-RF on the myobacterium tuberculosis data subset Lee2015
#' We use statistics of the SFS, of the minimal observable clade sizes and 
#' of the pairwise Hamming distances
#' as summary statistics for the ABC, using a random-forest based approach
#' We test the best-fitting model, Beta, vs. exponential decline

load("mtb_data_call90.RData")

i <- which(names(data_main)=="Lee2015")


rep1 <- 1
#' For setting random seeds, only partially applied
load("seeds_TB.RData")
seed1 <- seeds1[23*4+i+rep1]
set.seed(seed1)


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

#load("TB_data_ABC_data_fasta.RData")

#' number of cores for parallel computation
mc1 <- 7#20


gran1 <- c(log(0.5),log(5000))
gran2 <- c(log(0.5),log(20000))
f1 <- function(x){switch(x,gran1,gran2)}
growthrange0 <- lapply(c(1,2,1,1,2,1,2,2,1,1,1,2,2,1,1,1,2,2,2,1,1,1,1,1,2,2,2,2,2),f1)
betarange0 <- c(1,2)
diracrange0 <- c(0,1)
n0 <- 125000



#' Prior setting for exponential decline, uses 10K sims to 
#' approximate expected total tree length for Watterson's estimator
prior_obs_s_expd <- function(n_ind,nsimul=100,
                             ranges = c(-250,-200,-150,-100,-50,-25,-10),
                             s_obs = c(5,10,15,20,40,60,75)){
  watt_expdecl <- function(g){decline_start <- log(1000)/(-g)
  ms_com <- paste("-t 0.000001 -L -eG 0",g,"-eG",decline_start,"0")
  seq1 <- read.ms.output(ms(n_ind,10000,ms_com),
                         is.file = FALSE)
  return(mean(seq1$times[,2]))}
  theta <- sapply(ranges,function(g){if (g<0){th1 <- 2*s_obs/watt_expdecl(g)} else
  {th1 <- s_obs/sum((1:(n_ind-1))^(-1))}
    return(th1)})  
  sample_g <- sample(seq(along=ranges),nsimul,replace = TRUE)
  sample_s <- sample(seq(along=s_obs),nsimul,replace = TRUE)
  
  sample_theta <- rep(0,nsimul)
  for (i in 1:nsimul){sample_theta[i] <- theta[sample_s[i],sample_g[i]]}
  
  mod_vec <- rep(7,nsimul)
  coal_p <- ranges[sample_g]
  prior1 <- cbind(rep(n_ind,sum(nsimul)),mod_vec,coal_p,
                  s_obs[sample_s],sample_theta)
  colnames(prior1) <- c("n_ind","model","coal_param","s_obs","theta_watt")
  return(prior1)
}


#' Wrapper for simulation, including exponential decline (w. starting point trick)

sim_seq <- function(nsamp1,theta1,coal_param=0,model=1){
  KM <- (model==1 & coal_param==0) | (model==2 & length(coal_param)==1 & coal_param==2) | (model==6)
  KM <- KM | (model == 7 & coal_param==0)
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
  if (model==7 & !KM){
    decline_start <- log(1000)/(-coal_param)
    ms_com <- paste("-t",theta1,"-eG 0",coal_param,"-eG",decline_start,"0")
    raw_seq <- ms(nsam=nsamp1,opts=ms_com)
    seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])
  }
  
  if (sum(seq1)==0){return(NA)} 
  return(seq1)  
}

#' Diversity statistics
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

#' Get prior
setwd("../general_scripts/")
prior1 <- prior_obs_s_cont(data_main[[i]]$n_ind,models=c(2),nsimul=c(0,n0,0,0,0,0),
                      ranges = list(growthrange0[[i]],
                                    betarange0,diracrange0,NULL,NULL,0),
                      s_obs = rep(data_main[[i]]$n_mut,2),mc1 = mc1)

prior1 <- rbind(prior1,prior_obs_s_expd(data_main[[i]]$n_ind,nsimul = n0,
                                        s_obs = rep(data_main[[i]]$n_mut,2)))

prior1[,"theta_watt"] <- sapply(prior1[,"theta_watt"],log_smear)
setwd("../MTB_MMC_repo/")



clu1 <- makeForkCluster(nnodes = mc1)

sims1 <- parApply(clu1,prior1,1,function(x){
  divfun_ohaf(sim_seq(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],model = x[2]))})
stopCluster(clu1)

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

param1 <- coal_p[model1==7]
sumstat2 <- data.frame(param1,sumstat1[model1==7,])
growth_rf <- regAbcrf(param1~.,
                      data=sumstat2,
                      paral = TRUE)
growthfit_rf <- predict(growth_rf,target1,
                        sumstat2,paral=TRUE,
                        quantiles = quant_ret)


param1 <- coal_p[model1==2]
sumstat2 <- data.frame(param1,sumstat1[model1==2,])  
beta_rf <- regAbcrf(param1~.,
                    data=sumstat2,
                    paral = TRUE)
betafit_rf <- predict(beta_rf,target1,sumstat2,paral=TRUE,
                      quantiles = quant_ret)


model1 <- as.factor(model1)
sumstat1a <- data.frame(model1,t(sims1)) 

good_cols <- lapply(sumstat1a[,-1],
                    function(v){sum(tapply(v,sumstat1a[,1],function(x){range(x)[1]==range(x)[2]}))>1})
good_cols <- !(c(FALSE,unname(unlist(good_cols))))


rf_ms <- abcrf(model1~.,data=sumstat1a[,good_cols],paral = TRUE)

  
pred_ms <- predict(rf_ms,target1[good_cols[-1]],sumstat1a[,good_cols],
                   paral=TRUE)


asmmc_error <- 
  sum(rf_ms$model.rf$confusion.matrix["7",
                "2"])/sum(rf_ms$model.rf$confusion.matrix["7",c("7","2")])

modsel1 <- levels(pred_ms$allocation)[as.numeric(pred_ms$allocation)]
levels_rest <- levels(pred_ms$allocation)[-as.numeric(pred_ms$allocation)]
modsel2 <- levels_rest[which.max(pred_ms$vote[-as.numeric(pred_ms$allocation)])]
if (modsel1==6){postmed <- 0;post_q <- rep(NA,length(quant_ret))} else {
  fittedparam <- switch(modsel1,"1"=growthfit_rf,"2"=betafit_rf,"3"=diracfit_rf)
  postmed <- fittedparam$med; post_q <- fittedparam$quantiles
}
if (modsel2==6){postmed2 <- 0} else {
  fitted2param <- switch(modsel2,"7"=growthfit_rf,"2"=betafit_rf)
  postmed2 <- fitted2param$med 
}

out1 <- as.data.frame(as.list(c(dataset=names(data_main)[i],oob=round(rf_ms$prior.err,3),
                           asmmcoob=round(asmmc_error,3),
                      modelsel=modsel1,postp=round(pred_ms$post.prob,3),postmed = postmed, 
                      qu=post_q,modsel2=modsel2,postmed2=postmed2)))


write.table(out1,file=paste0("res/abcres",rep1,"_",names(data_main)[i],"expd.txt"),quote = FALSE,
            row.names = FALSE)




