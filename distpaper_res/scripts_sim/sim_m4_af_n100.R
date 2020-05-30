setwd("../")



set.seed(53373) #Seed 11 + 1, error with seed 53373: 
#Error in unserialize(node$con) : error reading from connection
#Calls: parApply ... FUN -> recvData -> recvData.SOCK0node -> unserialize
#Execution halted



#' Source scripts

source("../general_scripts/elength_lambdaxi.R")
source("../general_scripts/betaxicoal_sim.R")
source("../general_scripts/divstats.R")
source("construct_prior.R")


divfun_af <- function(seq1){
  if (is.matrix(seq1)){
    out1 <- c(quant_hm_oc(seq1),mean_sd_oc(seq1),
              hammfun(seq1),phylolength(seq1),r2fun(seq1),
              f_nucdiv_S(spectrum01(seq1)),allele_freqs(seq1))}
  else {out1 <- rep(NA,30)}
  names(out1) <-c("hm(O)",paste("O: qu",seq(.1,.9,.2)),"mean(O)","sd(O)",
                  paste("Ham: qu",seq(.1,.9,.2)),
                  paste("Phy: qu",seq(.1,.9,.2)),
                  paste("r2: qu",seq(.1,.9,.2)),
                  "Nucl. div.","S",
                  paste("AF: qu",seq(.1,.9,.2)))
  return(out1)}


library(parallel)

#' Argument list
nsamp <- 100
alpha_range <- seq(1,1.9,.1)
s_obs1 <- c(15,20,30,40,60,75)
#' Number of clusters
mc1 <- 7#16
#' We produce an additional replication of simulations 
#' Number of simulations per model class
nsim <- 175000



sim_seq <- function(nsamp1,theta1,coal_param=0,model=1){
  KM <- (model==1 & coal_param==0) | (model==2 & length(coal_param)==1 & coal_param==2) | (model==6)
  KM <- KM | (model==4 & length(coal_param)==1 & coal_param==2)
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
  if (model==4 & !KM){
    if (length(coal_param)==1){coal_param <- c(2-coal_param[1],coal_param[1])}
     seq1 <- betaxi_seq_sim(nsamp1,coal_param[1],coal_param[2],theta1)
  }
  if (sum(seq1)==0){return(NA)}
  return(seq1)
}

#' For time issues: Compute all possible rate functions.
#' Warning: Overwrites original rate function
#' \alpha\in seq(1,1.9,.1) is used here and n \leq 100
#' We make a lookup function


#beta4xi_rates_orig <- beta4xi_rates

#' Look up list: First coord. n, second alpha
#lookup_xirates <- vector("list",nsamp)

aux1 <- function(n){out1 <- lapply(alpha_range,function(a){
                                                    beta4xi_rates(n,2-a,a)})
                    names(out1) <- alpha_range
                    return(out1)
                    }
lookup_xirates <- mclapply(2:100,aux1,mc.cores = mc1)

#Rates as lookup function
beta4xi_rates <- function(n,a,b){if (a+b==2 & b %in% alpha_range){
                                 alph1 <- as.character(2-a)
                                 return(lookup_xirates[[n-1]][[alph1]])} else {
                                 return(NA)   
                                 }
                                 }

for (i in 1:2){

prior1 <- prior_obs_s(nsamp,models=c(4),nsimul=c(0,0,0,nsim,0,0),
              ranges = list(NULL,NULL,NULL,alpha_range,NULL,0),
              s_obs = s_obs1)
setwd("../distpaper_res/")

clu1 <- makeForkCluster(nnodes = mc1)

sims1 <- parApply(clu1,prior1,1,function(x){
  divfun_af(sim_seq(nsamp1 = x[1],theta1 = x[5],coal_param = x[3],model = x[2]))})

stopCluster(clu1)

save(prior1,sims1,file=paste0("sims_rep",i,"/sim_m4_af_n",nsamp,".RData"))
}
