#' Use on sims with stats O, Ham, AF+, r^2, Phy
#' Uses slightly different stat sets for ABC than abcinput_stand_stats.R
#' For abbreviations of stats, see Methods in our preprint

seed1 <- 46611 #Seed 14, use the same for all such ABCs

#' Rows that are kept for each set of statistics used for model selection
o_coarse <- 1:8
ham_coarse <- 9:13
phy_coarse <- 14:18
r2_coarse <-  19:23
s_pi <- c(24,25)
af_coarse <- 26:30
rows_abc <- list(c(o_coarse,ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse),
                 c(ham_coarse,phy_coarse),
                 c(o_coarse,r2_coarse),
                 c(s_pi,af_coarse),
                 c(o_coarse,s_pi,af_coarse),
                 c(ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse),
                 c(o_coarse,r2_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,phy_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,phy_coarse,r2_coarse),
                 c(ham_coarse,phy_coarse,s_pi,af_coarse))



#' Input vector
stats1 <- vector("list",10)

#' Names of sets of stats
names(stats1) <- c("FULL","Ham, Phy","O, r2","AF+","O, AF+",
                   "FULL - O","O, r2, AF+","FULL - r2","FULL - AF+","Ham,Phy,AF+")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",10)
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)
