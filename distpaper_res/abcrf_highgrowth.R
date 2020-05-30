#' Perform ABC only for "large" growth rates m16 
#' i.e. extract sims from both runs w. growth parameters > 20
#' and run abc (once)
#' sim_m16_af_n100.R needs to be run first

source("../general_scripts/misc_scripts.R") #' for debiased variable importance
library(abcrf)

seed1 <- 46611 #Seed 14, use the same for all such ABCs
set.seed(seed1)

prior_comb <- NULL
sims_comb <- NULL

for (k in 1:2){
load(paste0("sims_rep",k,"/sim_m16_af_n100.RData"))
prior_comb <- rbind(prior_comb,prior1)
sims_comb <- cbind(sims_comb,sims1)
}


high_g_rows <- which(prior_comb[,2]==1 & prior1[,3]>20)
no_g_rows <- sample(which(prior_comb[,2]==6),length(high_g_rows))
keep_rows <- c(high_g_rows,no_g_rows)


prior_comb <- prior_comb[keep_rows,]
sims_comb <- sims_comb[,keep_rows]

model1 <- prior_comb[,2]


#' Kick out any missing data (there should be none...occasionally there might be overflow)
bad_cols <- apply(sims_comb,2,function(v){any(is.na(v))})
bad_cols <- bad_cols 

sims_comb <- sims_comb[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
model1 <- as.factor(model1)

stats1 <- data.frame(model1,t(sims_comb))




#Exclude vars w/o variation
posvar_rvs <- 1+which(apply(stats1[,-1],2,sd)>0)
stay_cols <- c(1,posvar_rvs)
#' Use corrected variable importance  
rf1 <- abcrf_imp(model1~.,data=stats1[,stay_cols],ntree=500,paral = TRUE,ncores = 7, 
             lda=FALSE)
res_conf <- print(rf1)
res_imp <- variableImpPlot(rf1)



