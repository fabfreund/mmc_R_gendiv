#Comparing any two sets

oob_cf <- function(str1,str2){
load(paste0("abc_res/modsel_rep1_",str1,".RData"))
res_conf1 <- res_conf
load(paste0("abc_res/modsel_rep2_",str1,".RData"))
res_conf2 <- res_conf

oob1 <- sapply(res_conf1,
               function(l){sapply(l,
                                  function(m){mean(m[,ncol(m)])})})

oob2 <- sapply(res_conf2,
               function(l){sapply(l,
                                  function(m){mean(m[,ncol(m)])})})

oob_compare1 <-rbind(set1_rep1=apply(oob1,2,mean),set1_rep2=apply(oob2,2,mean))

load(paste0("abc_res/modsel_rep1_",str2,".RData"))
res_conf1 <- res_conf
load(paste0("abc_res/modsel_rep2_",str2,".RData"))
res_conf2 <- res_conf

oob1 <- sapply(res_conf1,
               function(l){sapply(l,
                                  function(m){mean(m[,ncol(m)])})})

oob2 <- sapply(res_conf2,
               function(l){sapply(l,
                                  function(m){mean(m[,ncol(m)])})})


oob_compare2 <-rbind(set2_rep1=apply(oob1,2,mean),
                     set2_rep2=apply(oob2,2,mean))
return(list(oob_compare1,oob_compare2))}


str1 <-"m12_smear2_n100"
str2 <- "m12_af_TB"
oob_cf(str1,str2)



