#' Use Rscript with 4+ arguments.  1: sims name (w/o sims_ and .RData), contains 
#' ALL sims needed (see next arguments) 
#' 2: model to test 3:,4:,... models to build forest off
args <- commandArgs(TRUE)

data_path <- paste0("sims_rep1/sim_",args[1],".RData")
#' mod1: this model's simulations will get sorted via a RF built on models mods_rf
mod1 <- as.numeric(args[2])
mods_rf <- as.numeric(args[-c(1,2)])

library(abcrf)

load(data_path)

#' Safety: Sims should not have missing data...
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})
bad_cols <- bad_cols 


sims1 <- sims1[,!(bad_cols)]
model1 <- as.factor(prior1[!(bad_cols),2])

#' sims for building decision RF
prior_m_rf <- which(model1 %in% mods_rf)
sims_m_rf <- sims1[,prior_m_rf]
model_m_rf <- droplevels(model1[prior_m_rf])
stats_m_rf <- data.frame(model_m_rf,t(sims_m_rf)) 

#' sims to allocate w. RF
prior_mtest <- which(model1==mod1)
sims_mtest <- sims1[,prior_mtest]
model_mtest <- droplevels(model1[prior_mtest])
stats_mtest <- data.frame(model_mtest,t(sims_mtest))


#' Kill stats that have no variation 
posvar_rvs <- 1+which(apply(stats_m_rf[,-1],2,sd)>0)
stay_cols <- c(1,posvar_rvs)

#' Build RF for model selection
rf_m_rf <- abcrf(model_m_rf~.,data=stats_m_rf[,stay_cols],ntree=1000,paral = TRUE,lda=FALSE)
#' For all sims of model mod1, allocate it to the other models via RF rf_m_rf
pred_m_rf <- predict(rf_m_rf,stats_mtest[,posvar_rvs],stats_m_rf[,stay_cols],paral = TRUE)
rm(rf_m_rf)

#' save model allocation result
save(pred_m_rf,file=paste0("mismod_res/pred_m",mod1,"vsrest_",args[1],".RData"))



