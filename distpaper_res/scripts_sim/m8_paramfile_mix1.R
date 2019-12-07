load("../seeds/seeds.RData")
set.seed(seeds1[20+3]) #Seed 20+i for mig scenario i (as in the filename)
#' Argument list
nsamp <- 100
expg_range <- c(0,0.5,1,2.5,4,7,10,25,50,75,100,500,1000)
s_obs1 <- c(75,75) #' Assumed observed mutations, always make two entries if
                   #' only one s_obs is used (we use the command <sample> to sample)

#' How subpopulations are subsampled (we assume equal-sized subpopulations, 
#' but not necessarily equal sampling)
popstruct <- list(c(10,90),c(20,80),c(30,70),c(60,40),c(50,50))

#' Range of migration rates
mig1 <- seq(0.1,2.1,0.5)


#' Number of simulations per model class per replication
nsim <- 175000/(length(mig1)*length(popstruct))

#' Number of clusters
mc1 <- 16

#' resulting name
simname <- "sim_m8_mix1.RData"
