args <- commandArgs(TRUE)

library(abcrf)
source("abcrf2.R")

for (k in 1:2){
load(paste0("../sims_rep",k,"/sim_",args[1],".RData"))


model1 <- prior1[,2]


stats1 <- vector("list",9)
res_conf <- vector("list",9)

names(stats1) <- c("full","noLDphylo","AF+","noOC","kato_rep","AF+OC","OC","OCHam","Ham")
names(res_conf) <- c("full","noLDphylo","AF+","noOC","kato_rep","AF+OC","OC","OCHam","Ham")


nreps <- 10
res_imp <- vector("list",nreps)


bad_cols <- apply(sims1,2,function(v){any(is.na(v))})
bad_cols <- bad_cols #| (prior1[,4]<=10)

sims1 <- sims1[,!(bad_cols)]
model1 <- model1[!(bad_cols)]
model1 <- as.factor(model1)

stats1[[1]] <- data.frame(model1,t(sims1)) #full
stats1[[2]] <- data.frame(model1,t(sims1[-(14:23),])) #All but LD + phylo branches
stats1[[3]] <- data.frame(model1,t(sims1[24:30,])) #Only SFS+
stats1[[4]] <- data.frame(model1,t(sims1[-(1:8),])) #all but OC
stats1[[5]] <- data.frame(model1,t(sims1[-c((1:8),24),])) #kato rep
stats1[[6]] <- data.frame(model1,t(sims1[-(9:23),])) #AF+OC
#stats1[[7]] <- data.frame(model1,t(sims1[c(1,24,25),])) # pi S hmOC
stats1[[7]] <- data.frame(model1,t(sims1[1:8,])) #OC 
#stats1[[9]] <- data.frame(model1,t(sims1[c(1:8,24,25),])) # pi S OC
stats1[[8]] <- data.frame(model1,t(sims1[1:13,])) # OC Hamming
stats1[[9]] <- data.frame(model1,t(sims1[9:13,])) # Hamming


for (i in 1:9){
res_conf[[i]] <- vector("list",nreps)  
print(i)
for (j in 1:nreps){
#' Model selection

#Exclude vars w/o variation
posvar_rvs <- 1+which(apply(stats1[[i]][,-1],2,sd)>0)
stay_cols <- c(1,posvar_rvs)
  
rf1 <- abcrf_imp(model1~.,data=stats1[[i]][,stay_cols],ntree=500,paral = TRUE,
             lda=FALSE)
res_conf[[c(i,j)]] <- print(rf1)
if (i==1){res_imp[[j]] <- variableImpPlot(rf1)}

}}



save(res_conf,res_imp,file=paste0("res_rep10_paper/modsel_rep",k,"_",args[1],"c.RData"))}
