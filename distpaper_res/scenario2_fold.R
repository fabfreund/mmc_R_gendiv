#' Name of simulated data (without sim_ and .RData) has to be both in sims_rep1 and sims_rep2)
dataname <- "m12_all_n100"
seed1 <- 60830 #Seed 13


sims_mod <- TRUE
names_mod <- TRUE
newname <- "m12_fold_n100"

sims_f <- function(m){
  #source fold function (sources too many funs...)
  source("../general_scripts/divstats.R")
  #fold the SFS
  temp1 <- apply(m[53:151,],2,fnl_sfs,scale1=FALSE)
  rownames(temp1) <- paste0("fS",1:nrow(temp1))
  # Kick out non-folded stats, add fSFS
  m <- rbind(m[-c(1:12,43:151),],temp1)}

#' Rows that are kept for each set of statistics used for model selection
ham_coarse <- c(1,3,5,7,9)
phy_coarse <- c(10,12,14,16,18)
r2_coarse <-  c(19,21,23,25,27)
s_pi <- c(28,29)
TajD <- 30
fsfs <- 31:80
  
rows_abc <- list(c(ham_coarse,phy_coarse,r2_coarse,s_pi,TajD,fsfs),
                 c(fsfs,s_pi,TajD),
                 r2_coarse,
                 ham_coarse,
                 phy_coarse,
                 c(r2_coarse,fsfs,s_pi,TajD),
                 c(ham_coarse,phy_coarse,fsfs,s_pi,TajD),
                 c(r2_coarse,ham_coarse,s_pi,TajD))



#' Input vector
stats1 <- vector("list",8)

#' Names of sets of stats
names(stats1) <- c("FULL","fSFS*","r2","Ham","Phy","r2, fSFS*",
                   "FULL - r2","r2, Ham,*")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",8)
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)
