#' Do model selection with ABC-RF, rejection method, 
#' multinomial logistic regression method and neural net regression method
#' Assess misclassification for a sample from the second set of simulations
#' Perform ABC based on the the first set of simulations
#' We need to have the simulations of models 1,2,3 (K+exp,Beta,Dirac) in
#' both sim folders

#ABC libraries
library(abc)
library(abcrf)
#ABC script
source("../general_scripts/misc_scripts.R")

library(parallel)


#' Number of clusters, number of sims taking for test purposes

mc1 <- 7
n_test <- 2500

#Take ntest sims with same params from second simulation run
load("sims_rep2/sim_m123_af_n100.RData")

bad_cols <- apply(sims1,2,function(v){any(is.na(v))})

sims1 <- sims1[,!(bad_cols)]
prior1 <- prior1[!(bad_cols),]

sample_rows <- c(sample(which(prior1[,"model"]==1),n_test),
                 sample(which(prior1[,"model"]==2),n_test),
                 sample(which(prior1[,"model"]==3),n_test))

prior_test <- prior1[sample_rows,]
sims_test <- sims1[,sample_rows] 


load("sims_rep1/sim_m123_af_n100.RData")

#Rename sims due to problems with original names (some special symbols...),
#maybe not necessary anymore...
rownames(sims_test)<- c("ochm",paste0("ocq",seq(1,9,2)),"ocmean","ocsd",paste0("hammq",seq(1,9,2)),
                             paste0("phyq",seq(1,9,2)),paste0("r2q",seq(1,9,2)),"pi","S",paste0("AFq",seq(1,9,2)))
                                           
rownames(sims1) <- c("ochm",paste0("ocq",seq(1,9,2)),"ocmean","ocsd",paste0("hammq",seq(1,9,2)),
                    paste0("phyq",seq(1,9,2)),paste0("r2q",seq(1,9,2)),"pi","S",paste0("AFq",seq(1,9,2)))



bad_cols <- apply(sims1,2,function(v){any(is.na(v))})

sims1 <- sims1[,!(bad_cols)]
prior1 <- prior1[!(bad_cols),]



zerovar_rvs <- which(apply(sims1,1,sd,na.rm=TRUE)==0)
sims1 <- sims1[-zerovar_rvs,]
sims_test <- sims_test[-zerovar_rvs,]

predmod_abc <- function(testa,priora,sima,tol1){
temp1 <- summary(postpr(target = testa,index =priora,sumstat = sima,tol = tol1,method = "neuralnet"))
res1 <- temp1$neuralnet$Prob
if (length(res1)<3){res1 <- c("1"=NA,"2"=NA,"3"=NA)}
names(res1) <- paste0("nnprob_",names(res1))
temp1  <- summary(postpr(target = testa,index =priora,sumstat = sima,tol = tol1,method = "mnlogistic"))
res2 <- temp1$mnlogistic$Prob
if (length(res2)<3){res2 <- c("1"=NA,"2"=NA,"3"=NA)}
names(res2) <- paste0("mnlprob_",names(res2))
res1 <- c(res1,res2)
res2 <- temp1$rejection$Prob
if (length(res2)<3){res2 <- temp1$Prob}
names(res2) <- paste0("rejprob_",names(res2))
res1 <- c(res1,res2)
return(res1)
}
  
clu1 <- makeForkCluster(nnodes = mc1)
test1 <- parApply(clu1,sims_test,2,predmod_abc,priora=prior1[,2],sima=t(sims1),tol1=0.005)
stopCluster(clu1)  

#load("ABC_allmeth.RData")

model1 <- as.factor(prior1[,2])
stats1 <- data.frame(model1,t(sims1))
#colnames(stats1) <- paste0("T",1:ncol(stats1))
#stats1 <- data.frame(model1=prior1[,2],t(sims1))
#rf1 <- abcrf_imp(model1~.,data=stats1,ntree=500,paral = TRUE,
#                 lda=FALSE)
rf1 <- abcrf(model1~.,data=stats1,ntree=500,paral = TRUE,
                                  lda=FALSE)



test_rf <- predict(object=rf1,obs=data.frame(t(sims_test)),
                   training=stats1,paral=TRUE)







miscl <- matrix(0,nrow=4,ncol=3,dimnames = list(c("rej","mnl","nn","rf"),
                                          c("MC EXP","MC BETA","MC Dirac")))

miscl_abc <- function(method,res1 = test1,priora =prior_test){
      miscl_counts <- rep(0,3)
      names(miscl_counts) <- c("1","2","3")
      for (i in 1:nrow(priora)){
      mod1 <- as.character(priora[i,2])
      if (method == "rej"){
      entry  <- switch(mod1,"1"=7,"2"=8,"3"=9)
      if (res1[entry,i]<=0.5){miscl_counts[mod1] <- miscl_counts[mod1] + 1} 
      }
      if (method == "mnl"){
        entry  <- switch(mod1,"1"=4,"2"=5,"3"=6)
        if (is.na(res1[entry,i])){entry  <- switch(mod1,"1"=7,"2"=8,"3"=9)}
        if (res1[entry,i]<=0.5){miscl_counts[mod1] <- miscl_counts[mod1] + 1} 
      }
      if (method == "nn"){
        entry  <- switch(mod1,"1"=1,"2"=2,"3"=3)
        if (is.na(res1[entry,i])){entry  <- switch(mod1,"1"=7,"2"=8,"3"=9)}
        if (res1[entry,i]<=0.5){miscl_counts[mod1] <- miscl_counts[mod1] + 1} 
      }
      }
      for (n1 in names(miscl_counts)){
        miscl_counts[n1] <- miscl_counts[n1]/table(priora[,2])[n1]}
      return(miscl_counts)}


miscl_abcrf <- function(res1 = test_rf,priora =prior_test){
      miscl_counts <- rep(0,3)
      names(miscl_counts) <- c("1","2","3")
      for (i in 1:nrow(priora)){
      mod1 <- as.character(priora[i,2])
      if (res1$allocation[i]!=mod1){
        miscl_counts[mod1] <- miscl_counts[mod1] + 1}} 
      for (n1 in names(miscl_counts)){
        miscl_counts[n1] <- miscl_counts[n1]/table(priora[,2])[n1]}
      return(miscl_counts)
      }
miscl[1,] <- miscl_abc("rej")
miscl[2,] <- miscl_abc("mnl")
miscl[3,] <- miscl_abc("nn")
miscl[4,] <- miscl_abcrf()

save(miscl,file="ABC_allmeth_m123_altrf.RData")
