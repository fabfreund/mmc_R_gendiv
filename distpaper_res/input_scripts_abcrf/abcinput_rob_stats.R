#' Use on sims with stats nO, nHam, nSFSl
#' For abbreviations of stats, see Methods in our preprint

seed1 <- 92210 #Seed 19

#' Rows that are kept for each set of statistics used for model selection
o_coarse <- 1:5
ham_coarse <- 6:10
SFSl <- 11:25
rows_abc <- list(o_coarse,
                 ham_coarse,
                 SFSl,
                 c(ham_coarse,SFSl),
                 c(o_coarse,ham_coarse),
                 c(o_coarse,SFSl),
                 c(o_coarse,ham_coarse,SFSl))
#' Input vector
stats1 <- vector("list",7)
#' Which stat set is used for importance scoring?
index_imp <- 7
#' Names of sets of stats
#names(stats1) <- c("FULL","nO","nHam","nSFSl",
#                   "nHam,nSFSl","nO,nSFSl","nO,nHam")
names(stats1) <- c("nO","nHam","nSFSl",
                   "nHam,nSFSl","nO,nHam","nO,nSFSl","FULL")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",7)
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)
