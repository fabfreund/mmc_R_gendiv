#' Reduce sims from model classes K+exp (=1) and Beta (=2) to
#' the sims from Dirac (=3). Sims are AF+, r^2, Phy, Ham and O (see preprint)
#' Glue all together 
#for (k in 1:2){
#load(paste0("sims_rep",k,"/sim_m3_af_n100.RData"))
#prior_temp <- prior1
#sims_temp <- sims1
#load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
#sims1 <- sims1[which(rownames(sims1) %in% rownames(sims_temp)),]
#prior1 <- rbind(prior1,prior_temp)
#sims1 <- cbind(sims1,sims_temp)
#save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
#}

#' same for BetaXi = m4
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m4_af_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
  sims1 <- sims1[which(rownames(sims1) %in% rownames(sims_temp)),]
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m124_af_n100.RData"))
}




#' Glue Dirac+exp = m5 with m123
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  load(paste0("sims_rep",k,"/sim_m5_af_n100.RData"))
  prior1 <- prior1[,which(colnames(prior1) %in% colnames(prior_temp)),]
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m1235_af_n100.RData"))
}
#' Glue Dirac+exp = m5 vers2 with m123
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  load(paste0("sims_rep",k,"/sim_m5_af_n100_largeg.RData"))
  prior1 <- prior1[,which(colnames(prior1) %in% colnames(prior_temp)),]
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m1235_af_n100_largeg.RData"))
}

#Glue m3 and m5
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m3_af_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  load(paste0("sims_rep",k,"/sim_m5_af_n100.RData"))
  prior1 <- prior1[,which(colnames(prior1) %in% colnames(prior_temp)),]
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m35_af_n100.RData"))
}

#Glue m3 and m5 vers 2
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m3_af_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  load(paste0("sims_rep",k,"/sim_m5_af_n100_largeg.RData"))
  prior1 <- prior1[,which(colnames(prior1) %in% colnames(prior_temp)),]
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m35_af_n100_largeg.RData"))
}

#' Glue m8 variants with m12 - all stats
for (set1 in c("p5p5","p9p1"))
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  if (set1=="p5p5"){adda <- ""} else {adda <- "b"}
  load(paste0("sims_rep",k,"/sim_m8_",set1,"_n100",adda,".RData"))
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m128_all",set1,"_n100.RData"))
}

#' Glue m8 variants with m12 
for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m8_mix1.RData"))
  prior_temp <- prior1
  sims_temp <- sims1
  load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
  prior1 <- rbind(prior1,prior_temp)
  sims1 <- cbind(sims1,sims_temp)
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m128_allmix1_n100.RData"))
}


#' Glue ALL m8 variants with m128_allmix1 and only keep set af... 
#' Define first which stats are at which row
o_coarse <- c(1,2,4,6,8,10:12)
ham_coarse <- seq(13,21,2)
phy_coarse <- seq(22,30,2)
r2_coarse <- seq(31,39,2)  
s_pi <- c(40,41)
af_coarse <- seq(44,52,2)
keep_stats <- c(o_coarse,ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse)

for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
  prior_temp <- prior1
  sims_temp <- sims1[keep_stats,]
  for (set1 in c("p5p5","p9p1")){
  if (set1=="p5p5"){adda <- "";newmod <- 81} else {adda <- "b";newmod <- 82}
  load(paste0("sims_rep",k,"/sim_m8_",set1,"_n100",adda,".RData"))
  prior1[,2] <- newmod
  sims1 <- sims1[keep_stats,]
  prior_temp <- rbind(prior1,prior_temp)
  sims_temp <- cbind(sims1,sims_temp)}
  for (set1 in c("mix1","mix2")){
    if (set1=="mix1"){newmod <- 83} else {newmod <- 84}
  load(paste0("sims_rep",k,"/sim_m8_",set1,".RData"))
  prior1[,2] <- newmod
  sims1 <- sims1[keep_stats,]
  prior_temp <- rbind(prior1,prior_temp)
  sims_temp <- cbind(sims1,sims_temp)}
  prior1 <- prior_temp
  sims1 <- sims_temp
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m12all8_af.RData"))
}

#' Also add reciprocal m1 and m8 sims for last approach, just for sims_rep1/

load("sims_rep2/sim_m12_all_n100.RData")
prior1[prior1[,2]==1,2] <- 11 
prior1[prior1[,2]==2,2] <- 12 
prior1 <- prior1[prior1[,2] %in% c(11,12),]
sims_temp <- sims1[keep_stats,prior1[,2] %in% c(11,12)]
prior_temp <- prior1
load("sims_rep1/sim_m12all8_af.RData")
prior1 <- rbind(prior1,prior_temp)
sims1 <- cbind(sims1,sims_temp)
save(prior1,sims1,file="sims_rep1/sim_m12all8m12b_af.RData")
#for (k in 1:2){
#  load(paste0("sims_rep",k,"/sim_m8_mix1.RData"))
#  prior_temp <- prior1
#  sims_temp <- sims1
#  load(paste0("sims_rep",k,"/sim_m12_all_n100.RData"))
#  prior1 <- rbind(prior1,prior_temp)
#  sims1 <- cbind(sims1,sims_temp)
#  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m128_allmix1_n100.RData"))
#}


#Extract simulation sets of pairs of model classes from m123
for (mod1 in c(1,2)){
for (k in 1:2){
load(paste0("sims_rep",k,"/sim_m123_af_n100.RData"))
keep1 <- (prior1[,"model"] %in% c(mod1,3))   
prior1 <- prior1[keep1,]
sims1 <- sims1[,keep1]
save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m",mod1,"3_af_n100.RData"))
}
}  


#Extract pairwise simulation sets of model classes from m124
#for (mod1 in c(1,2)){
  for (k in 1:2){
    load(paste0("sims_rep",k,"/sim_m124_af_n100.RData"))
    keep1 <- (prior1[,"model"] %in% c(2,4))   
    prior1 <- prior1[keep1,]
    sims1 <- sims1[,keep1]
    save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m24_af_n100.RData"))
  }

for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m124_af_n100.RData"))
  keep1 <- (prior1[,"model"] %in% c(1,4))   
  prior1 <- prior1[keep1,]
  sims1 <- sims1[,keep1]
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m14_af_n100.RData"))
}  


#Extract pairwise simulation sets of model classes from m128
for (set1 in c("p5p5","p9p1")){
  for (mod1 in 1:2){
  for (k in 1:2){
  load(paste0("sims_rep",k,"/sim_m128_all",set1,"_n100.RData"))
  keep1 <- (prior1[,"model"] %in% c(mod1,8))   
  prior1 <- prior1[keep1,]
  sims1 <- sims1[,keep1]
  save(prior1,sims1,file=paste0("sims_rep",k,"/sim_m",mod1,"8_",set1,"_all.RData"))
}}}