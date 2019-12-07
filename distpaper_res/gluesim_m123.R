#' Reduce sims from model classes K+exp (=1) and Beta (=2) to
#' the sims from Dirac (=3). Sims are AF+, r^2, Phy, Ham and O (see preprint)
for (k in 1:2){
load(paste0("sims_rep",k,"/sim_m3_af_n100.RData"))
prior_temp <- prior1
sims_temp <- sims1
load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
sims1 <- sims1[which(rownames(sims1) %in% rownames(sims_temp)),]
prior1 <- rbind(prior1,prior_temp)
sims1 <- cbind(sims1,sims_temp)
save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
}

#Extract simulation sets of pairs of model classes
#for (mod1 in c(1,2)){
#for (k in 1:2){
#load(paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
#keep1 <- (prior1[,"model"] %in% c(mod1,3))   
#prior1 <- prior1[keep1,]
#sims1 <- sims1[,keep1]
#save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m",mod1,"3_af_n100.RData"))
#}
#}  