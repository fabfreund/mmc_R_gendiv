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

#' Example (not run, and only works if simulations and ABC was run before...)
#' Important: the full ABC results are e.g. named as modsel_rep2_m12_sfsl_n100.RData, the function only needs the middle part  
#' str1 <-"m12_fold_n100" 
#' str2 <- "m12_sfsl_n100"
#' oob_cf(str1,str2)

#' Report mean misclassification matrices
#' input is the ABC main name (w/o modsel_rep1 or _rep2 and .RData)
extract_conf <- function(data1){
extract <- function(l){  
out1 <- l[[1]]
for (i in 2:length(l)){out1 <- out1 + l[[i]]}
out1 <- out1[,-ncol(out1)]
rowsums1 <- rowSums(out1)
for (i in 1:nrow(out1)){
  out1[i,] <- out1[i,]/rowsums1[i]}
return(out1)
}
load(paste0("abc_res/modsel_rep1_",data1,".RData"))
res_conf1 <- res_conf
load(paste0("abc_res/modsel_rep2_",data1,".RData"))
extract_set <- 
function(set1){(extract(res_conf1[[set1]])+extract(res_conf[[set1]]))/2}
out1 <- lapply(names(res_conf),extract_set)
names(out1) <- names(res_conf)
return(out1)
}
