set.seed(27589) #Seed 3
#' Source scripts


source("../general_scripts/lambdacoal_sim.R")
source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/ext_funs.R")
source("../general_scripts/divstats.R")
source("construct_prior.R")


#' Wrapper for the different simulation scripts
#' We use ms to simulate Kingman and related coalescents

library(phyclust) #R implementation of Hudson's ms
library(gap) #Read in output of ms
#library(parallel)

sim_seq <- function(nsamp1,theta1,coal_param=0,model=1){
  KM <- (model==1 & coal_param==0) | (model==2 & length(coal_param)==1 & coal_param==2) | (model==6)
  if (KM){raw_seq <- ms(nsam=nsamp1,opts=paste("-t",theta1))
  seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])}
  if (model==1 & !KM){
    raw_seq <- ms(nsam=nsamp1,opts=paste("-t",theta1,"-G",coal_param))
    seq1 <- t(read.ms.output(raw_seq,is.file = FALSE)$gametes[[1]])
  }
  if (model==2 & !KM){
    if (length(coal_param)==1){coal_param <- c(2-coal_param[1],coal_param[1])}
    if (all(coal_param == c(1,1))){seq1 <- bsz_seq_sim(nsamp1,theta1)} else {
      seq1 <- beta_seq_sim(nsamp1,coal_param[1],coal_param[2],theta1)}
  }
  if (model==3){
    seq1 <-  dirac_seq_sim(nsamp1,coal_param,theta1)
  }
  if (sum(seq1)==0){return(NA)}
  return(seq1)
}




times1 <- vector("list",5)
names(times1) <- c("Phy","Ham","O","LD","AF")
times2 <- times1;times3 <- times1 

for (rep1 in 1:10){
#' Do 1000 priors, switch to general_scripts to use the expected length C script
curdir <- getwd()
setwd("../general_scripts/")
prior1 <- prior_obs_s(100,models=c(1,2),nsimul=c(1000,1000,0,0,0,0),
                      ranges = list(c(0,0.5,1,2.5,4,7,10,25,50,75,100,500,1000),
                                    seq(1,1.9,0.1),NULL,NULL,NULL,0),
                      s_obs = c(15,20,30,40,60,75))

prior2 <- prior_obs_s(n_ind = 147,models = c(1,2),
                      nsimul=c(1000,1000,0,0,0,0),
                      ranges = list(c(0,0.5,1,2.5,4,7,10,25,50,75,100,500,1000),
                                    seq(1,1.9,0.1),NULL,NULL,NULL,0),
                      s_obs = c(454,454)
)

prior3<- prior_obs_s(n_ind = 250,models = c(1,2),
                      nsimul=c(1000,1000,0,0,0,0),
                      ranges = list(c(0,0.5,1,2.5,4,7,10,25,50,75,100,500,1000),
                                    seq(1,1.9,0.1),NULL,NULL,NULL,0),
                      s_obs = c(500,500)
)
setwd(curdir)

seq1 <- apply(prior1,1,function(x){sim_seq(nsamp1 = x[1],
                                           theta1 = x[5],
                                           coal_param = x[3],model = x[2])})
seqTB <- apply(prior2,1,function(x){sim_seq(nsamp1 = x[1],
                                            theta1 = x[5],
                                            coal_param = x[3],model = x[2])})
seq3 <- apply(prior3,1,function(x){sim_seq(nsamp1 = x[1],
                                           theta1 = x[5],
                                           coal_param = x[3],model = x[2])})

f_O <- function(seq2){if (is.matrix(seq2)){quant_hm_oc(seq2)}}
f_Ham <- function(seq2){if (is.matrix(seq2)){hammfun(seq2)}}
f_Phy <- function(seq2){if (is.matrix(seq2)){phylolength(seq2)}}
f_AF <- function(seq2){if (is.matrix(seq2)){c(f_nucdiv_S(spectrum01(seq2)),
                                              allele_freqs(seq2))}}
f_LD <- function(seq2){if (is.matrix(seq2)){r2fun(seq2)}}


times1[[1]] <- rbind(times1[[1]],system.time(lapply(seq1,f_Phy)))
times1[[2]] <- rbind(times1[[2]],system.time(lapply(seq1,f_Ham)))
times1[[3]] <- rbind(times1[[3]],system.time(lapply(seq1,f_O)))
times1[[4]] <- rbind(times1[[4]],system.time(lapply(seq1,f_LD)))
times1[[5]] <- rbind(times1[[5]],system.time(lapply(seq1,f_AF)))

times2[[1]] <- rbind(times2[[1]],system.time(lapply(seqTB,f_Phy)))
times2[[2]] <- rbind(times2[[2]],system.time(lapply(seqTB,f_Ham)))
times2[[3]] <- rbind(times2[[3]],system.time(lapply(seqTB,f_O)))
times2[[4]] <- rbind(times2[[4]],system.time(lapply(seqTB,f_LD)))
times2[[5]] <- rbind(times2[[5]],system.time(lapply(seqTB,f_AF)))

times3[[1]] <- rbind(times3[[1]],system.time(lapply(seq3,f_Phy)))
times3[[2]] <- rbind(times3[[2]],system.time(lapply(seq3,f_Ham)))
times3[[3]] <- rbind(times3[[3]],system.time(lapply(seq3,f_O)))
times3[[4]] <- rbind(times3[[4]],system.time(lapply(seq3,f_LD)))
times3[[5]] <- rbind(times3[[5]],system.time(lapply(seq3,f_AF)))
}

save(times1,times2,times3,file="times_divstats.RData")  

