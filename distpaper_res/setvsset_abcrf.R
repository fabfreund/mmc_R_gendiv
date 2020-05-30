#' Script to test how good the OOB errors match the real error 
#' from an independent data set
#' We test this for K+exp vs B vs Dirac three-fold comparison

#' Load first data set
data_path <- "sims_rep1/sim_m123_af_n100.RData"

library(abcrf)

load(data_path)

#' Safety: Sims should not have missing data...
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})
bad_cols <- bad_cols 
sims1 <- sims1[,!(bad_cols)]
model1 <- as.factor(prior1[!(bad_cols),2])

#' sims for building decision RF
sims_m_rf <- sims1
model_m_rf <- model1
stats_m_rf <- data.frame(model_m_rf,t(sims_m_rf)) 

#' sims to allocate w. RF
load("sims_rep2/sim_m123_af_n100.RData")

#' Safety: Sims should not have missing data...
bad_cols <- apply(sims1,2,function(v){any(is.na(v))})
bad_cols <- bad_cols 
sims1 <- sims1[,!(bad_cols)]
model1 <- as.factor(prior1[!(bad_cols),2])


sims_mtest <- sims1
model_mtest <- model1
stats_mtest <- data.frame(model_mtest,t(sims_mtest))

rm(prior1,sims1)

#' Kill stats that have no variation 
posvar_rvs <- 1+which(apply(stats_m_rf[,-1],2,sd)>0)
stay_cols <- c(1,posvar_rvs)

#' Build RF for model selection
#rf_m_rf <- abcrf(model_m_rf~.,data=stats_m_rf[,stay_cols],ntree=1000,paral = TRUE,
#                 ncores = 7, 
#                 lda=FALSE)
#' Build RF for model selection with unbiased variable importances
#' Will slightly alter errors compared to the command above
rf_m_impc <- abcrf_imp(model_m_rf~.,data=stats_m_rf[,stay_cols],ntree=1000,
                       ncores = 7,
                       paral = TRUE,lda=FALSE)

#' For all sims of the independent test set, 
#' allocate them to the model classes via decision RF built by the other sim set
#pred_m_rf <- predict(rf_m_rf,stats_mtest[,posvar_rvs],stats_m_rf[,stay_cols],paral = TRUE)
#te <- cbind(model_mtest,pred_m_rf$allocation)
#error_kexp <- sum(te[te[,1]==1,1]!=te[te[,1]==1,2])/nrow(te[te[,1]==1,])
#error_beta <- sum(te[te[,1]==2,1]!=te[te[,1]==2,2])/nrow(te[te[,1]==2,])
#error_dirac <- sum(te[te[,1]==3,1]!=te[te[,1]==3,2])/nrow(te[te[,1]==3,])


pred_m_imp <- predict(rf_m_impc,stats_mtest[,posvar_rvs],stats_m_rf[,stay_cols],paral = TRUE)
te <- cbind(model_mtest,pred_m_imp$allocation)
error_kexp_ic <- sum(te[te[,1]==1,1]!=te[te[,1]==1,2])/nrow(te[te[,1]==1,])
error_beta_ic <- sum(te[te[,1]==2,1]!=te[te[,1]==2,2])/nrow(te[te[,1]==2,])
error_dirac_ic <- sum(te[te[,1]==3,1]!=te[te[,1]==3,2])/nrow(te[te[,1]==3,])

#rm(rf_m_rf)
#conf_m <- rf_m_rf$model.rf$confusion.matrix
conf_m_impc <- rf_m_impc$model.rf$confusion.matrix

#' save model allocation result
#save(conf_m,error_kexp,error_beta,error_dirac,pred_m_rf,
#     file="mismod_res/simvssim_m123_af_n100.RData")
save(error_kexp_ic,error_beta_ic,error_dirac_ic,pred_m_imp,
     conf_m_impc,
     file="mismod_res/simvssim_m123_af_n100_impc.RData")


