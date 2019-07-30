#' Name of simulated data (without sim_ and .RData) has to be both in sims_rep1 and sims_rep2)
dataname <- "scenario1"
seed1 <- 50654 #Seed 6

#' Rows that are kept for each set of statistics used for model selection
rows_abc <- list(1:30,-(14:23),24:30,-(1:8),-c((1:8),24),-(9:23),1:8,1:13,9:13)



#' Input vector
stats1 <- vector("list",9)

#' Names of sets of stats
names(stats1) <- c("FULL","O, Ham, AF+","AF+","FULL - O","FULL - (O,pi)","O, AF+","O","O, Ham","Ham")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",9)
names(res_conf) <- c("FULL","O, Ham, AF+","AF+","FULL - O","FULL - (O,pi)","O, AF+","O","O, Ham","Ham")
res_imp <- vector("list",nreps)
