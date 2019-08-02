#' Name of simulated data (without sim_ and .RData) has to be both in sims_rep1 and sims_rep2)
dataname <- "m12_af_TB"
seed1 <- 57285 #Seed 7

#' Rows that are kept for each set of statistics used for model selection
rows_abc <- list(1:30,-(14:23),-(9:23),-(19:23),-(14:18),-(9:13),-(1:8))



#' Input vector
stats1 <- vector("list",7)

#' Names of sets of stats
names(stats1) <- c("FULL","O, Ham, AF+","O, AF+","FULL - r2","FULL - PHY","FULL - Ham","FULL - O")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",7)
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)
