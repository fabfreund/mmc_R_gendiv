#' Use on sims with O, Ham, r^2, AF+
#' In paper, used for sims without singleton mutations
#' For abbreviations of stats, see Methods in our preprint

seed1 <- 79309 #Seed 18

#' Rows that are kept for each set of statistics used for model selection
o_coarse <- 1:8
ham_coarse <- 9:13
r2_coarse <- 14:18
s_pi <- c(19,20)
af_coarse <- 21:25
rows_abc <- list(c(o_coarse,ham_coarse,r2_coarse,s_pi,af_coarse),
                 c(ham_coarse),
                 c(o_coarse,r2_coarse),
                 c(s_pi,af_coarse),
                 c(ham_coarse,r2_coarse,s_pi,af_coarse),
                 c(o_coarse,r2_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,r2_coarse),
                 c(ham_coarse,s_pi,af_coarse))



#' Input vector
stats1 <- vector("list",9)

#' Names of sets of stats
names(stats1) <- c("FULL","Ham","O, r2","AF+",
                   "FULL - O","FULL - Ham","FULL - r2","FULL - AF+","Ham, AF+")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",9)
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)
