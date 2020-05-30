load("../seeds/seeds.RData")
#set.seed(seeds1[20+2]) #Seed 20+i for mig scenario i (as in the filename)
set.seed(seeds1[20+3]) #Seed 20+i for mig scenario i (as in the filename)
#' Argument list
nsamp <- 100
expg_range <- c(0,0.5,1,2.5,4,7,10,25,50,75,100,500,1000)
s_obs1 <- c(75,75) #' Assumed observed mutations, always make two entries if
                   #' only one s_obs is used (we use the command <sample> to sample)

#' How subpopulations are subsampled (we assume equal-sized subpopulations, 
#' but not necessarily equal sampling)
popstruct <- list(c(90,10))

#' Range of migration rates
mig1 <- c(0.5)


#' Number of simulations per model class per replication
nsim <- 175000

#' Number of clusters
mc1 <- 7#16

#' resulting name
#simname <- "sim_m8_p9p1_n100.RData"
simname <- "sim_m8_p9p1_n100b.RData"