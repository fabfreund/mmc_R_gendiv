#' Expected total length E(Ln) of Lambda and Beta-Xi n-coalescent.
#' Computing recursively via M\"ohle '06, Adv. in Appl. Probab., Eq. 2.3


#Input: nmax sample size for which to compute E(Ln), ratefun Function of n that
#gives the transition rates of the block counting pocess of the Lambda-n-coalescent. 
#If s > 0 it returns the
#Watterson estimator 2*s/E(Ln) instead.
eln_f_lambda <- function(nmax,ratefun,s=0){
eln <- rep(0,nmax)
rates_l <- sapply(2:nmax,ratefun) #list of Rate vectors for all n used
eln[2] <- 2/(sum(rates_l[[c(1,2)]])) #E(L_2)=2*E(T2)
for (i in 3:nmax){
    rates1 <- rates_l[[i-1]] #entry i-1 contains rate for i blocks 
    #Compute transition probabilities for a certain jump size
    #A jump from i to j blocks means coalescing i-j+1 blocks
    p_i <- rep(0,i-1)
    for (j in 1:(i-1)){	
      p_i[j] <- rates1[i-j+1] #rate from i to j
    }
    #trans probs = rates/totalrate
    p_i <- p_i/(sum(rates1))
    #Recursion
    eln[i] <- i/(sum(rates1)) + sum(eln[1:(i-1)]*p_i)}
  if (s>0){eln[nmax] <- 2*s/eln[nmax]}
  return(eln[nmax])  
}



#' Compute E(Ln) for the BetaXi n-coalescent
#' Input: nmax sample size for which to compute E(Ln), ratefun beta4xi_rates
#' from script betaxicoal_sim.R as a function of n, e.g. function(n){beta4xi_rates(n,0.5,1.5)}
#' for BetaXi with params a=0.5,b=1.5.
#If s > 0 it returns the
#Watterson estimator 2*s/E(Ln) instead.

eln_f_xi <- function(nmax,ratefun,s=0){
  eln <- rep(0,nmax)
  rates_l <- lapply(2:nmax,ratefun) #Rate vectors for all n used
  eln[2] <- 2/(sum(rates_l[[c(1,2)]])) #E(L_2)=2*E(Exp(0.25))
  for (i in 3:nmax){
    rates1 <- rates_l[[i-1]] 
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
  if (s>0){eln[nmax] <- 2*s/eln[nmax]}
  return(eln[nmax])
}
