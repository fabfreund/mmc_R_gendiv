#' Name of simulated data (without sim_ and .RData) has to be both in sims_rep1 and sims_rep2)
dataname <- "m12_all_n100"
seed1 <- 25966 #Seed 12


sims_mod <- TRUE
names_mod <- TRUE
newname <- "m12_sfsl_n100"

sims_f <- function(m){m <- rbind(m,"S15+"=colSums(m[-(1:66),]))}

#' Rows that are kept for each set of statistics used for model selection
o_coarse <- c(1,2,4,6,8,10,11,12)
ham_coarse <- c(13,15,17,19,21)
phy_coarse <- c(22,24,26,28,30)
r2_coarse <-  c(31,33,35,37,39)
s_pi <- c(40,41)
HD <- c(42,43)
sfsl <- c(53:66,152)
jere <- c(53,152)
  
rows_abc <- list(c(o_coarse,ham_coarse,phy_coarse,r2_coarse,s_pi,HD,sfsl),
                 c(s_pi,HD,sfsl),
                 sfsl,
                 jere,
                 c(o_coarse,r2_coarse,jere),
                 c(o_coarse,r2_coarse,sfsl),
                 c(o_coarse,r2_coarse,s_pi),
                 c(r2_coarse,sfsl),
                 c(o_coarse,ham_coarse,s_pi,jere)
                 )



#' Input vector
stats1 <- vector("list",9)

#' Names of sets of stats
names(stats1) <- c("FULL","SFSl*","SFSl","Sing.t","O, r2, Sing.t","O, r2, SFSl",
                   "O, r2, S, pi","r2,SFSl","O, Ham, Sing.t+")


#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and varuiable importances
res_conf <- vector("list",9)
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)
