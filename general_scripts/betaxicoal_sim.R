#(Fourfold)Beta-Xi coalescent simulator
#Jump chain partition, jump chain block counting, sequences 

#Auxiliary functions
#Multinomial coefficient: Argument tuple k_vect of integers. Output multinomial coefficient
#sum(k_vect) choose entries of k_vect.
multinom_hb <- function(k_vect){k <- sum(k_vect)
temp1 <- 1
if (length(k_vect)>1){ #if FALSE, then we have k scalar with k=n in the main sc$
  for (i in 1:(length(k_vect)-1)){
    temp1 <- temp1*choose(sum(k_vect[i:length(k_vect)]),k_vect[i])
  }}

return(temp1)
}


#Compute the potential merger counts and their total rates 
#following Blath et al "The site-frequency spectrum associated with \Xi-coalescents", 
#Eq. 15 (corrected), for the BetaXi (their rescaled version)
#Note that there is a typo in Eq. 15 of Blath et al, the power of 4 is -(k+l), see Eq 14 therein 
#Compute Eq. 15 for a certain combination of sum and number of merged blocks (this is all it depends on)
#For arguments, see Eq. 15
eq15 <- function(beta_p1,beta_p2,no_blocks,sum_merged,number_mergers){
  a <- beta_p1; b<- beta_p2; m <- no_blocks; k <- sum_merged; r <- number_mergers  
  sum_max <- min(m-k,4-r)
  M_fact <- c(4,12,24,24) #all possible values for the falling factorial in Eq. 14
  res1 <- 0
  for (l in 0:sum_max){
    res1 <- res1 + choose(m-k,l)*M_fact[r+l]*4^(-k-l)*beta(a+k+l-2,b+m-k-l)  
  }
  return(res1/beta(a,b))
}

#Transition rates BetaXi for every distinct set of sizes for block sets to be merged to
#block sets. As an example, we compute  we have one rate for all mergers of two sets of 
#two blocks each (merger class {2,2}). Input n: sample size, a,b: Beta-Xi parameters 
beta4xi_rates <- function(n,a,b){
  #build the possible merger classes
  sets <- NULL
  sets <- as.list(2:n)
  temp1 <- sets
  for (m in 2:4){
    temp3 <- NULL
    for (k in temp1){
      if (sum(k)<= n-2){  
        temp2 <- as.list(2:min(n-sum(k),min(k)))
        temp2 <- lapply(temp2,append,k,after=0)
        temp3 <- c(temp3,temp2)
        sets <- c(sets,temp2)
        temp1 <- temp3}}
  }
  #compute rate for each merger class by summing Eq. 14
  rates1 <- rep(0,length(sets))
  for (s in seq(along=sets)){
    k_sum <- sum(sets[[s]])
    r  <- length(sets[[s]]) 
    k <- sets[[s]]
    #Multinomial coeff. counts multiple same-sized picks, this is correcting that
    multiply_corr <- prod(sapply(unname(table(k)),factorial))
    if (k_sum!=n) {k <- c(k,n-k_sum)}
    rates1[s] <- multinom_hb(k)*eq15(a,b,n,k_sum,r)/multiply_corr 
  }
  return(list(sets,rates1))
}

#Simulate the jump chain of a Beta-Xi n-coalescent
#
#Arguments: n=sample size, 
#ratef=transition rates of the Beta-Xi n-coalescent, as computed above,
#as a funtion of n
#
#Output: Matrix, row i gives states of the jump chain at jump i-1 of the Beta-Xi n-coalescent
#Row i shows integers that indicate the partition of [n] at this jump
#Blocks are defined by who has the same integer

Xi_jc_sim <- function(n,ratef){config <- matrix(-(1:n),nrow=1) #At jump 0, everybody is in its own block
                              jumpcount <- 1 #Actually jumps+1
                              mergecount <- 0
                              repeat{
                              blocks <- unique(config[jumpcount,]) #present blocks before jump   
                              n1 <- length(blocks) #number of blocks
                              k_configs <- ratef(n1)[[1]]
                              rates1 <- ratef(n1)[[2]] 
                              transprob <- rates1/sum(rates1) #Transition probabilities
                              k <- unlist(sample(k_configs,1,prob = transprob)) #which block config is chosen?
                              newblock <- vector("list",length(k))
                              temp_block <- sample(blocks) #Makes a random permutation of the blocks
                              #Which blocks are merged: 
                              #Just take the first k1, then the next k2,... of temp_block
                              l1 <- 0 #How many blocks are taken already for mergers
                              for (j in 1:length(k)){newblock[[j]] <- temp_block[(l1+1):(l1+k[[j]])]
                                                     l1 <- l1 + k[[j]]}
                              temp <- config[jumpcount,] #The block configuration before jump
                              #merge blocks:
                              jumpcount <- jumpcount + 1
                              for (j in 1:length(k)){mergecount <- mergecount+1
                              temp[temp %in% newblock[[j]]] <- mergecount} 
                              config <- rbind(config,temp)
                              if (length(unique(config[jumpcount,]))==1){break}
                              }
                              return(unname(config))}

#Simulate sequences according to a given jump chain
#
#Arguments: conf1 is output from Xi_jc_sim, theta scaled mutation rate,
#ratef is the transition rate function of Beta-Xi, as above
#
#Output: individuals in rows, mutations in columns
Xi_seq_sim <- function(conf1,theta,ratef){
                        totratef <- function(n){sum(ratef(n)[[2]])} #Xi-rate functions return a list, 2nd entry 
                        n <- ncol(conf1) 
                        seq1 <- NULL
                        theta <- theta/2 #Mutations with half rate s.t. E(\pi_{12})=\theta 
                        bc_jc <- apply(conf1,1,function(x){length(unique(x))}) #Get block counts
                        for (i in 1:(nrow(conf1)-1)){
                          blocks <- unique(conf1[i,]) #Blocks present at jump
                          n1 <- length(blocks) #number of blocks
                          tot_rate <- totratef(n1) #total rate with n1 blocks
                          #mut1 <- rgeom(n1,tot_rate/(tot_rate+theta)) #number mut's is geometric
                          wt <- rexp(1,tot_rate)
                          mut1 <- rpois(n1,theta*wt)
                          for (j in 1:n1){
                            if (mut1[j]>0){
                            temp <- rep(0,n)
                            temp[conf1[i,]==blocks[j]] <- 1
                            temp <- matrix(temp,ncol=mut1[j],nrow=n)
                            seq1 <- cbind(seq1,temp)}
                          }
                        }                 
                        return(seq1)}

#Wrapper script for generating DNA sequences under the Beta-Xi n-coalescents
#Arguments as above


betaxi_seq_sim <- function(n,a,b,theta){ratef1 <- function(m){beta4xi_rates(m,a,b)}
                                      Xi_seq_sim(Xi_jc_sim(n,ratef1),theta,ratef1)}


