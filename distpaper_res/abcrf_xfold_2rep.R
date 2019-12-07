#' Run via Rscript. Provide 2 or 3 arguments
#' Argument 1: ABC input script (see subfolder input_scripts_abcrf/)
#' Argument 2: sims name (without "sim_" and ".RData", from sims_rep1 and sims_rep2)
#' Argument 3: Only when argument 1 is abcinput_all.R, then sets which option
#' in abcinput_all.R is used (see there)

args1 <- commandArgs(TRUE)

lda1 <- FALSE # Should linear discriminants be included in the ABC as summ. stats.?
change_name <- FALSE # FALSE if input argument should determine output name 
                     # (see last line) 

source("../general_scripts/misc_scripts.R") #' for debiased variable importance
library(abcrf)

seed1 <- 46611 #Seed 14, use the same for all such ABCs
set.seed(seed1)
#' For which set of stats do we record variable importances (usually set 1 = FULL)
#' sets of stats are defined in the input script (argument 1)
index_imp <- 1

for (k in 1:2){
dataname <- args1[2]
load(paste0("sims_rep",k,"/sim_",dataname,".RData"))
source(paste0("input_scripts_abcrf/",args1[1]))

if (lda1){lda_params <- vector("list",nreps)}

model1 <- prior1[,2]


#' Kick out any missing data (there should be none...occasionally there might be overflow)
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})
bad_cols <- bad_cols 

sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
model1 <- as.factor(model1)

#' Extract sets of stats 
for (i in 1:length(stats1)){stats1[[i]] <- data.frame(model1,t(sims1[rows_abc[[i]],]))}


for (i in 1:length(stats1)){
res_conf[[i]] <- vector("list",nreps)  
print(i)
for (j in 1:nreps){
#' Model selection (nreps times) for set i of statistics

#Exclude vars w/o variation
posvar_rvs <- 1+which(apply(stats1[[i]][,-1],2,sd)>0)
stay_cols <- c(1,posvar_rvs)
#' Use corrected variable importance  
rf1 <- abcrf_imp(model1~.,data=stats1[[i]][,stay_cols],ntree=500,paral = TRUE, 
             lda=lda1)
res_conf[[c(i,j)]] <- print(rf1)
if (i==index_imp){res_imp[[j]] <- variableImpPlot(rf1)
if (lda1){lda_params[[j]] <- rf1$model.lda}
}

}}

if (change_name){dataname <- paste0(dataname,args1[3])}


save(res_conf,res_imp,file=paste0("abc_res/modsel_rep",k,"_",dataname,".RData"))
if (lda1){save(lda_params,file=paste0("abc_res/lda_rep",k,"_",dataname,".RData"))}
}
