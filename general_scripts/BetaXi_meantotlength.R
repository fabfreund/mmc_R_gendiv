#' Recursion for expected total length of the Beta-Xi coalescent
#' POLISHING NOT FINISHED 

#' Main problem: Computing BetaXi rates is costly, and we do this for every
#' Watterson estimate and simulation. To bypass this, we now compute a table
#' for the rates and then use the table instead of calling the function anew
library(parallel)
source("betaxicoal_sim.R")
mc1 <- 3
n_start <- as.numeric(args[2])
alpha_range <- seq(1,1.9,0.1)
rates_l <- vector("list",length(alpha_range))
for (i in seq(along=alpha_range)){
  b1 <- alpha_range[i] #alpha
  a1 <- 2 - b1 #2-alpha  
  rates_l[[i]] <- mclapply(2:n_start,beta4xi_rates,a=a1,b=b1,mc.cores = mc1,
                           mc.preschedule = FALSE)
}

#' Recursion for E(Ln), nmax is n until which to recursively compute
#' Temporarily: alpha_pos is one of the entries of alpha_range
eln_f_betaxi <- function(nmax,alpha_pos){
  eln <- rep(0,nmax)
  eln[2] <- 2/(sum(rates_l[[c(alpha_pos,1,2)]])) #E(L_2)=2*E(Exp(0.25))
  for (i in 3:nmax){
    rates1 <- rates_l[[c(alpha_pos,i-1)]] 
    #Compute transition probabilities for a certain jump size
    #e.g.: A collision of k1,k2 blocks leads to the state n-sum(k)+2 (2 new 
    #blocks...)
    #states1 is reduction of blocks
    states1 <- sapply(rates1[[1]],sum)-sapply(rates1[[1]],length)  
    p_i <- rep(0,i-1)
    for (j in 1:(i-1)){	
      p_i[j] <- sum(rates1[[2]][which(states1 == i-j)])
    }
    #trans probs = rates/totalrate
    p_i <- p_i/(sum(rates1[[2]]))
    eln[i] <- i/(sum(rates1[[2]])) + sum(eln[1:(i-1)]*p_i)}  
  return(eln[nmax])
}