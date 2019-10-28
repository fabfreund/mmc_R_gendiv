#' read in argument as name of input file

args1 <- commandArgs(TRUE)

lda1 <- FALSE
source("../general_scripts/misc_scripts.R") #' for debiased variable importance
library(abcrf)

seed1 <- 46611 #Seed 14, use the same for all such ABCs
set.seed(seed1)


for (k in 1:2){
dataname <- args1[2]
load(paste0("sims_rep",k,"/sim_",dataname,".RData"))
source(args1[1])

if (lda1){lda_params <- vector("list",nreps)}

model1 <- prior1[,2]


#' Kick out any missing data (there should be none...occasionally there might be overflow)
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})
bad_cols <- bad_cols #| (prior1[,4]<=10)

sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
model1 <- as.factor(model1)

for (i in 1:length(stats1)){stats1[[i]] <- data.frame(model1,t(sims1[rows_abc[[i]],]))}


for (i in 1:length(stats1)){
res_conf[[i]] <- vector("list",nreps)  
print(i)
for (j in 1:nreps){
#' Model selection

#Exclude vars w/o variation
posvar_rvs <- 1+which(apply(stats1[[i]][,-1],2,sd)>0)
stay_cols <- c(1,posvar_rvs)
#' Use corrected variable importance  
rf1 <- abcrf_imp(model1~.,data=stats1[[i]][,stay_cols],ntree=500,paral = TRUE, 
             lda=lda1)
res_conf[[c(i,j)]] <- print(rf1)
if (i==1){res_imp[[j]] <- variableImpPlot(rf1)
if (lda1){lda_params[[j]] <- rf1$model.lda}
}

}}

if (change_name){dataname <- paste0(dataname,args1[3])}


save(res_conf,res_imp,file=paste0("abc_res/modsel_rep",k,"_",dataname,".RData"))
if (lda1){save(lda_params,file=paste0("abc_res/lda_rep",k,"_",dataname,".RData"))}
}
