#' Script to simulate SNPs from a Beta coalescent genealogy (infinite sites model)
#' Uses a setting close to Eldholm 2015 (but without serial sampling)


#' Seed setting
load("seeds_TB.RData")
seed1 <- seeds1[12]
set.seed(seed1)

#' Step I: load Eldholm data set, report properties
#' We will use the frequencies to draw the bases for the mock data
library(pegas)
seq_eld <- read.dna("data_fasta/Eldholm2015.fasta",
                    format = "fasta")
bfreq_eld <- base.freq(seq_eld[-1,]) #Omit ancestral sequence, record ACTG freqs

#' Function to recode 1-0 SNP matrices to mock DNA (bases drawn from real ACTG 
#' distribution) 
mock_dna <- function(v){nucbase <- sample(x = names(bfreq_eld),
                                          size = 2,prob = bfreq_eld)
                        v1 <- rep("z",length(v))
                        v1[v==1] <- nucbase[1]
                        v1[v==0] <- nucbase[2]
                        return(v1)}
#' load simulation scripts: Watterson estimator computation and simulation script

source("../general_scripts/lambdacoal_sim.R")
source("../general_scripts/elength_lambdaxi.R")

#' Simulation parameters
#alpha <- c(1,1.5,1.25,1.9)
alpha <- c(0.1,0.25,0.5,0.75,1.75)
nsamp <- 250
sobs <- 500
nsim <- 10

#' Produce simulations 
for (a in alpha){
sims1 <- vector("list",length = nsim)
ratef1 <- function(n){beta_bc_rates(n,a = 2-a,b=a)}
watttheta <- eln_f_lambda(nmax = nsamp,ratefun = ratef1,s=sobs)
for (i in 1:nsim){
sim_temp <- beta_seq_sim(n = nsamp,a = 2-a,b = a,theta = watttheta)
sim_transf <- apply(sim_temp,2,mock_dna) #Transform SNP matrices
write.dna(sim_transf,file = paste0("sim_fasta/sim_eldlike_alpha",
                                     a,"_rep",i,".fasta"),nbcol = -1,
          format = "fasta",colsep = "")
sims1[[i]] <- sim_transf  
}

save(sims1,file=paste0("sims_eldlike_beta_",a,".RData"))    
}
