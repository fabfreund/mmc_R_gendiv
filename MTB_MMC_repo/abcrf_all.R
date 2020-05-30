#' Script to perform analyses 1-4,6-8 (see Sheet 1 in Sup_table_2.xlsx, 
#' supplementary to the biorxiv preprint, version 3)
#' Run via Rscript with three integer arguments
#' 1st arg is data set number from mtb_data_call*.RData, 2nd arg is seed set, 
#' 3rd argument is number of cores used for analysis
#' Prerequisite: Run for the data set specified by the 1st argument
#' sim_all_sets.R

args <- commandArgs(TRUE)
i <- as.integer(args[1])
rep1 <- as.integer(args[2])
mc1 <- as.integer(args[3])

load(paste0("sims/fullsims_rep",rep1,"_dataset",i,".RData"))

#' Set seed as if main analysis would be 
#.Random.seed <- sim_seed

#' R packages 
library(abcrf)
library(parallel)
library(psych)


#' Scripts from this repo
source("../general_scripts/misc_scripts.R")
source("../general_scripts/divstats.R")



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




run_abc <- function(prior1,sims1,target1){

model1 <- prior1[,2]
coal_p <- prior1[,3]

#' For safety: Remove simulations with NA (due to overflow)
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})

sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
coal_p <- coal_p[!(bad_cols)]

sumstat1 <- t(sims1)
rownames(sumstat1) <- NULL
 
#quantiles to be returned for posterior distribution
quant_ret <- seq(0,1,0.025)
#' Parameter estimation Kingman+exp
param1 <- coal_p[model1==1]
sumstat2 <- data.frame(param1,sumstat1[model1==1,])
growth_rf <- regAbcrf(param1~.,
                      data=sumstat2,
                      paral = TRUE,ncores = mc1)
growthfit_rf <- predict(growth_rf,target1,
                        sumstat2,paral=TRUE,
                        ncores=mc1,
                        quantiles = quant_ret)

#' Parameter estimation Beta
param1 <- coal_p[model1==2]
sumstat2 <- data.frame(param1,sumstat1[model1==2,])  
beta_rf <- regAbcrf(param1~.,
                    data=sumstat2,
                    paral = TRUE,ncores = mc1)
betafit_rf <- predict(beta_rf,target1,sumstat2,paral=TRUE,
                      ncores=mc1,
                      quantiles = quant_ret)

#' Parameter estimation Dirac
param1 <- coal_p[model1==3]
sumstat2 <- data.frame(param1,sumstat1[model1==3,])   
dirac_rf <- regAbcrf(param1~.,
                     data=sumstat2,
                     paral = TRUE,ncores = mc1)
diracfit_rf <- predict(dirac_rf,target1,sumstat2,
                       paral=TRUE,ncores=mc1,
                       quantiles = quant_ret)



model1 <- as.factor(model1)
sumstat1a <- data.frame(model1,t(sims1)) 
#' Remove summary statistics with no variation within a model class
good_cols <- lapply(sumstat1a[,-1],
                    function(v){sum(tapply(v,sumstat1a[,1],function(x){range(x)[1]==range(x)[2]}))>=1}) #added that single no var classes already cause problems...
good_cols <- !(c(FALSE,unname(unlist(good_cols))))

#' ABC model selection
rf_ms <- abcrf(model1~.,data=sumstat1a[,good_cols],paral = TRUE,
               ncores=mc1)
pred_ms <- predict(rf_ms,target1[good_cols[-1]],sumstat1a[,good_cols],
                   paral=TRUE,ncores=mc1,
                   paral.predict=TRUE,
                   ncores.predict=mc1)

#Build error output
asmmc_error <- 
  mean(c(sum(rf_ms$model.rf$confusion.matrix["1",c("2","3")])/sum(rf_ms$model.rf$confusion.matrix["1",c("1","2","3","6")]),
         sum(rf_ms$model.rf$confusion.matrix["6",c("2","3")])/sum(rf_ms$model.rf$confusion.matrix["6",c("1","2","3","6")])))

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

out1 <- as.data.frame(as.list(c(dataset=names(data_main)[i],oob=round(rf_ms$prior.err,3),
                           asmmcoob=round(asmmc_error,3),
                      modelsel=modsel1,postp=round(pred_ms$post.prob,3),postmed = postmed, 
                      qu=post_q,modsel2=modsel2,postmed2=postmed2)))
return(list(out1,growthfit_rf,betafit_rf,diracfit_rf,pred_ms))
}





#' 1) Load data 90 or 70
#' 2) Extract stats of the data set
#' 3) Perform all ABC, main analysis is seeded as before 


#' Load the TB data 90 
load("mtb_data_call90.RData")
target1 <- as.data.frame(t(divfun_ohaf(data_main[[i]]$seq_0_1)))
.Random.seed <- sim_seed
#' KM+exp unif prior, Beta unif on [1,2], Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_90[[1]],prior_90[[2]],
                      prior_90[[3]],prior_90[[4]]),
                sims1 = cbind(sims_90[[1]],sims_90[[2]],
                              sims_90[[3]],sims_90[[4]]),
                target1 = target1)

addon <- "unig_beta1_redstat_het90"
write.table(res1[[1]],file=paste0("res/abcres",rep1,"_",names(data_main)[i],
                             "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/res",rep1,names(data_main)[i],addon,".RData"))

#' KM+exp log prior, Beta unif on [1,2], Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_90[[5]],prior_90[[2]],
                               prior_90[[3]],prior_90[[4]]),
                sims1 = cbind(sims_90[[5]],sims_90[[2]],
                              sims_90[[3]],sims_90[[4]]),
                target1 = target1)

addon <- "logg_beta1_redstat_het90"
write.table(res1[[1]],file=paste0("res/abcres",rep1,"_",names(data_main)[i],
                                  "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/res",rep1,names(data_main)[i],addon,".RData"))


#' KM+exp log prior, Beta unif on [0,2], Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_90[[5]],prior_90[[6]],
                               prior_90[[3]],prior_90[[4]]),
                sims1 = cbind(sims_90[[5]],sims_90[[6]],
                              sims_90[[3]],sims_90[[4]]),
                target1 = target1)

addon <- "logg_beta0_redstat_het90"
write.table(res1[[1]],file=paste0("res/abcres",rep1,"_",names(data_main)[i],
                                  "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/res",rep1,names(data_main)[i],addon,".RData"))

#' KM+exp log prior, Beta unif on [0,2], 1% replaced by BSZ, 
#' Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_90[[5]],prior_90[[7]],
                               prior_90[[3]],prior_90[[4]]),
                sims1 = cbind(sims_90[[5]],sims_90[[7]],
                              sims_90[[3]],sims_90[[4]]),
                target1 = target1)

addon <- "logg_beta0_redstat_het90"
write.table(res1[[1]],file=paste0("res_bsz/abcres_bsz",rep1,
                                  "_",names(data_main)[i],
                                  "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata_bsz/res",rep1,names(data_main)[i],addon,".RData"))

#' Load the TB data 75 
load("mtb_data_call75.RData")
target1 <- as.data.frame(t(divfun_ohaf(data_main[[i]]$seq_0_1)))
.Random.seed <- sim_seed
#' KM+exp unif prior, Beta unif on [1,2], Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_75[[1]],prior_75[[2]],
                               prior_75[[3]],prior_75[[4]]),
                sims1 = cbind(sims_75[[1]],sims_75[[2]],
                              sims_75[[3]],sims_75[[4]]),
                target1 = target1)

addon <- "unig_beta1_redstat_het75"
write.table(res1[[1]],file=paste0("res/abcres",rep1,"_",names(data_main)[i],
                                  "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/res",rep1,names(data_main)[i],addon,".RData"))

#' KM+exp log prior, Beta unif on [1,2], Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_75[[5]],prior_75[[2]],
                               prior_75[[3]],prior_75[[4]]),
                sims1 = cbind(sims_75[[5]],sims_75[[2]],
                              sims_75[[3]],sims_75[[4]]),
                target1 = target1)

addon <- "logg_beta1_redstat_het75"
write.table(res1[[1]],file=paste0("res/abcres",rep1,"_",names(data_main)[i],
                                  "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/res",rep1,names(data_main)[i],addon,".RData"))


#' KM+exp log prior, Beta unif on [0,2], Dirac, KM
res1 <- run_abc(prior1 = rbind(prior_75[[5]],prior_75[[6]],
                               prior_75[[3]],prior_75[[4]]),
                sims1 = cbind(sims_75[[5]],sims_75[[6]],
                              sims_75[[3]],sims_75[[4]]),
                target1 = target1)

addon <- "logg_beta0_redstat_het75"
write.table(res1[[1]],file=paste0("res/abcres",rep1,"_",names(data_main)[i],
                                  "_",addon,".txt"),quote = FALSE,
            row.names = FALSE)
growthfit_rf <- res1[[2]]
betafit_rf <- res1[[3]]
diracfit_rf <- res1[[4]]
pred_ms <- res1[[5]]

save(growthfit_rf,diracfit_rf,betafit_rf,pred_ms,
     file=paste0("resdata/res",rep1,names(data_main)[i],addon,".RData"))
