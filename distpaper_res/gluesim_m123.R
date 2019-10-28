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


for (mod1 in c(1,2)){
for (k in 1:2){
load(paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
keep1 <- (prior1[,"model"] %in% c(mod1,3))   
prior1 <- prior1[keep1,]
sims1 <- sims1[,keep1]
save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m",mod1,"3_af_n100.RData"))
}
}  