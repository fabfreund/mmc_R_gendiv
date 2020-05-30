#' Simulation of time-serial sampling in Lambda-coalescents and K+exp
#' Different sampling times, allow time changes
#' Coal_times: sampling times in coalescent time units 

#' The code below is commented out, but fully functional. 
#' We provide a revised version below which has better functionality
#' Additional packages
#library(extraDistr)
#Lambda_difft_coal <- function(coal_times,ratef1,kmexp=FALSE,rho=1){
#  rho <- rho/2 #To match the 4N0 scaling of Hudson's ms
#  ratef <- function(m){if (m>1){out1 <- ratef1(m)} else {out1 <- 0}
#    return(out1)}
#  coal_times <- unname(coal_times)
#  coal_counts <- as.vector(table(coal_times))
#  coal_times2 <- as.numeric(labels(table(coal_times))[[1]])
#  n <- length(coal_times) # How many sequences are there in total
#  config <- list(-(1:coal_counts[1])) #At jump 0,time 0 ind. are in their own blocks
#  time1 <- 1
#  i <- 1
#  jt <- NULL 
#  n_untnow <- coal_counts[1] 
#  tottime <- 0
#  repeat{
#    blocks <- unique(config[[i]]) #present blocks before jump   
#    n1 <- length(blocks)  #number of blocks
#    if ((time1 == length(coal_times2)) & (n1 == 1)) {break}
#    #breaks if we have sampled all ind. AND only a single block remains 
#    wt1 <- rexp(1,sum(ratef(n1))) #Waiting time for next coalescence
#    #So far only implemented for Kingman: 
#    #add growth following Tavare,Griffiths p. 404
#    #sampling theory..., 1994
#    #We use rgompertz from extraDistr, here we need to set
#    #a <- b*eta, b is our b
#    if (kmexp){
#      gomp_a <- choose(n1,2)*exp(rho*tottime)
#      gomp_b <- rho
#      wt1 <- rgompertz(1,gomp_a,gomp_b)}
#    time_cond <- FALSE
#    if (time1  < length(coal_times2)){
#      if ((n1==1) | (wt1 > coal_times2[time1+1]-tottime)) {time_cond <- TRUE}}
#    if (time_cond){
#      time1 <- time1+1;jt[i] <- coal_times2[time1]-tottime; 
#      tottime <- tottime + jt[i];i <- i+1 
#      config[[i]] <- c(config[[i-1]],-(n_untnow+(1:coal_counts[time1])))
#      n_untnow <- n_untnow+coal_counts[time1]} else {
#        jt[i] <- wt1; tottime <- tottime + jt[i] 
#        transprob <- ratef(n1)/sum(ratef(n1)) #Transition probabilities
#        k <- sample(1:n1,1,prob = transprob) #How many blocks are merged?
#        newblock <- sample(blocks,k) #Which blocks are merged
#        temp <- config[[i]] #From the block configuration before jump
#        i <- i+1 
#        temp[temp %in% newblock] <- i #merge blocks
#        config[[i]] <- temp
#      }}
#  out2 <- list(jc=unname(config),coalt = jt)
#  return(out2)}

#' REVISED VERSION of the serial coalescent, adds a general time change function
#' G(t) (more precisely its inverse function) 
#' as e.g. in Eq. 22 "Cannings models, population size changes and multiple merger coalescents"
#' (Freund, JOMB, 2020, see also Eq. 10).
#' Time change is provided as argument Ginv which is an one-dimensional function 
#' and the inverse of the time change function from Freund 2020
#' You need to determine G defined as in the mentioned paper yourself
#' Default is no time-change G(t)=t
#' If you use it for Kingman and compare it to ms, use half growth rate rho/2 
#' if you want to match the 4N0 scaling of Hudson's ms

Lambda_difft_coal_tc <- function(coal_times,ratef1,Ginv=function(t){t}){
  ratef <- function(m){if (m>1){out1 <- ratef1(m)} else {out1 <- 0}
    return(out1)}
  coal_times <- unname(coal_times)
  coal_counts <- as.vector(table(coal_times))
  coal_times2 <- as.numeric(labels(table(coal_times))[[1]])
  n <- length(coal_times) # How many sequences are there in total
  config <- list(-(1:coal_counts[1])) #At jump 0,time 0 ind. are in their own blocks
  time1 <- 1
  i <- 1
  jt <- NULL 
  n_untnow <- coal_counts[1] 
  tottime <- 0
  repeat{
    blocks <- unique(config[[i]]) #present blocks before jump   
    n1 <- length(blocks)  #number of blocks
    if ((time1 == length(coal_times2)) & (n1 == 1)) {break}
    #breaks if we have sampled all ind. AND only a single block remains 
    wt1 <- rexp(1,sum(ratef(n1))) #Waiting time for next coalescence
    time_cond <- FALSE
    if (time1  < length(coal_times2)){
      if ((n1==1) | (wt1 > coal_times2[time1+1]-tottime)) {time_cond <- TRUE}}
    if (time_cond){
      time1 <- time1+1;jt[i] <- coal_times2[time1]-tottime; 
      tottime <- tottime + jt[i];i <- i+1 
      config[[i]] <- c(config[[i-1]],-(n_untnow+(1:coal_counts[time1])))
      n_untnow <- n_untnow+coal_counts[time1]} else {
        jt[i] <- wt1; tottime <- tottime + jt[i] 
        transprob <- ratef(n1)/sum(ratef(n1)) #Transition probabilities
        k <- sample(1:n1,1,prob = transprob) #How many blocks are merged?
        newblock <- sample(blocks,k) #Which blocks are merged
        temp <- config[[i]] #From the block configuration before jump
        i <- i+1 
        temp[temp %in% newblock] <- i #merge blocks
        config[[i]] <- temp
      }}
  jt2 <- rep(-1,length(jt))
  jt2[1] <- Ginv(jt[1]) 
  if (length(jt)>1){
  for (i in seq(along=jt)[-1]){jt2[i] <- Ginv(sum(jt[1:i])) - Ginv(sum(jt[1:(i-1)]))}}
  out2 <- list(jc=unname(config),coalt = jt2)
  return(out2)}


#' Simulate sequences according to a given coalescent
#' conf1 is output from Lambda_difft_coal
#' seq1 (output) shows individuals in rows, mutations in columns

Lambda_seq_sim_difft<- function(coal1,theta){
  conf1 <- coal1$jc
  coalt <- coal1$coalt
  n <- max(sapply(conf1,length))
  seq1 <- NULL
  theta <- theta/2 #Mutations with half rate s.t. E(\pi_{12})=\theta 
  for (i in 1:(length(conf1)-1)){
    blocks <- unique(conf1[[i]]) #Blocks present at jump
    n1 <- length(blocks) #number of blocks #total rate with n1 blocks
    mut1 <- rpois(n1,theta*coalt[i]) #number mut's given branch length is Poissonian
    for (j in 1:n1){
      if (mut1[j]>0){
        temp <- rep(0,length(conf1[[i]]))
        temp[conf1[[i]]==blocks[j]] <- 1
        temp <- c(temp,rep(0,n-length(temp)))
        temp <- matrix(temp,ncol=mut1[j],nrow=n)
        seq1 <- cbind(seq1,temp)}
    }
  }                 
  return(seq1)}

#' Function to get total length from output of  Lambda_difft_coal
#' Just multiplicate number of blocks present (=number diff. numbers 
#' in chain) w. waiting time. Ignore last state of chain (=MRCA) 

length_serialcoal <- function(coal1){
  nblocks <- sapply(coal1$jc,
                    function(x){length(unique(x))})
  nblocks <- nblocks[-length(nblocks)]
  wait_times <- coal1$coalt
  return(sum(nblocks*wait_times))}

#' same for height
height_serialcoal <- function(coal1){
  wait_times <- coal1$coalt
  return(sum(wait_times))}

#' use 1K sims to approximate expected total length 
#' (just 1K due to computational cost)
est_elength_serial <- function(ct,rate1,kme,rho,mc1){
  h1 <- function(x){
    length_serialcoal(Lambda_difft_coal(ct,rate1,kme,rho))}
  ret1 <- unlist(mclapply(1:1000,h1,mc.cores = mc1))  
  return(mean(ret1))
}

est_elength_serial_tc <- function(ct,rate1,Ginv,mc1,mcsched=TRUE){
  h1 <- function(x){
    length_serialcoal(Lambda_difft_coal_tc(ct,rate1,Ginv))}
  ret1 <- unlist(mclapply(1:1000,h1,mc.cores = mc1,
                          mc.preschedule = mcsched))  
  return(mean(ret1))
}

#' Same for height
est_eheight_serial <- function(ct,rate1,kme,rho,mc1){
  h1 <- function(x){
    height_serialcoal(Lambda_difft_coal(ct,rate1,kme,rho))}
  ret1 <- unlist(mclapply(1:1000,h1,mc.cores = mc1))  
  return(mean(ret1))
}

est_eheight_serial_tc <- function(ct,rate1,Ginv,mc1){
  h1 <- function(x){
    height_serialcoal(Lambda_difft_coal_tc(ct,rate1,Ginv))}
  ret1 <- unlist(mclapply(1:1000,h1,mc.cores = mc1))  
  return(mean(ret1))
}



