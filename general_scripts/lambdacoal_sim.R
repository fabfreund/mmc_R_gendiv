#Lambda-Coalescent simulator

#Simulate the jump chain of a Lambda n-coalescent
#
#Arguments: n=sample size, 
#ratef=transition probabilities of the jump chain of the Lambda n-coalescent 
#as a funtion of n
#
#Output: Matrix, row i gives states of the jump chain at jump i-1 of the Lambda n-coalescent
#Row i shows integers that indicate the partition of [n] at this jump
#Blocks are defined by who has the same integer

Lambda_jc_sim <- function(n,ratef){config <- matrix(-(1:n),nrow=1) #At jump 0, everybody is in its own block
i <- 1
repeat{
  blocks <- unique(config[i,]) #present blocks before jump   
  n1 <- length(blocks) #number of blocks
  transprob <- ratef(n1) #Transition probabilities
  k <- sample(1:n1,1,prob = transprob) #How many blocks are merged?
  newblock <- sample(blocks,k) #Which blocks are merged
  temp <- config[i,] #From the block configuration before jump
  i <- i+1 
  temp[temp %in% newblock] <- i #merge blocks
  config <- rbind(config,temp)
  if (length(unique(config[i,]))==1){break}
}
return(unname(config))}

#Simulate sequences according to a given jump chain
#
#Arguments: conf1 is output from Lambda_jc_sim, theta scaled mutation rate,
#totratef: the total rate for the next jump as a funtion of n
#
#Output: individuals in rows, mutations in columns
Lambda_seq_sim <- function(conf1,theta,totratef){
  n <- ncol(conf1) 
  seq1 <- NULL
  theta <- theta/2 #Mutations with half rate s.t. E(pw. diff. of 2 seqs)=\theta 
  bc_jc <- apply(conf1,1,function(x){length(unique(x))}) #Get block counts
  for (i in 1:(nrow(conf1)-1)){
    blocks <- unique(conf1[i,]) #Blocks present at jump
    n1 <- length(blocks) #number of blocks
    tot_rate <- totratef(n1) #total rate with n1 blocks
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



#Wrapper scripts for generating DNA sequences under Dirac, Beta, Bolthausen-Sznitman-n
#coalescents
#
#First, transition probabilities and total rates
#and rate that exactly k blocks are merged, 
#For transiion probabilities and rates: first entry is 0 since 
#rates/probabilities are available for 2,3,... blocks


#For Bolthausen-Sznitman n-coalescent. n=sample size
#Transition probabilities 
bsz_bc_transprob <- function(n){x <- n/((n-1)*(2:n)*(1:(n-1)))
return(c(0,x))}
#Total rate function
bsz_bc_rate <- function(n){n-1}

#Rates block counting of Beta(a,b)-n-coalescent. n=sample size
beta_bc_rates <-function(n,a,b){x<-rep(0,n)
if (a==0 & b==2){x[2] <- choose(n,2)} else {
  for (k in 2:n) {x[k]<-choose(n,k)*beta(a+k-2,n-k+b)/beta(a,b)}}
return(x)}
#Transition probabilities 
beta_bc_transprob <- function(n,a,b){x<-beta_bc_rates(n,a,b)
return(x/sum(x))}
#Total rate function
beta_bc_rate<- function(n,a,b){sum(beta_bc_rates(n,a,b))}

#Rates block counting of Dirac-n-coalescent, Lambda point mass in p. n=sample size 
dirac_bc_rates <-function(n,p){x<-rep(0,n)
for (k in 2:n) {x[k]<-choose(n,k)*p^(k-2)*(1-p)^(n-k)}
return(x)}
#Transition probabilities
dirac_bc_transprob <- function(n,p){x <- dirac_bc_rates(n,p)
return(x/sum(x))}
#Total rate
dirac_bc_rate<- function(n,p){sum(dirac_bc_rates(n,p))}

#Wrapper simulate n sequences Bolthausen-Sznitman n-coalescent
bsz_seq_sim <- function(n,theta){Lambda_seq_sim(Lambda_jc_sim(n,bsz_bc_transprob),theta,
                                                bsz_bc_rate)}   
#Wrapper simulate n sequences Beta(a,b)-n-coalescent
beta_seq_sim <- function(n,a,b,theta){beta_trans <- function(n){beta_bc_transprob(n,a,b)}
beta_rate <- function(n){beta_bc_rate(n,a,b)}
Lambda_seq_sim(Lambda_jc_sim(n,beta_trans),theta,
               beta_rate)}

#Wrapper simulate n sequences Dirac(p)-n-coalescent
dirac_seq_sim <- function(n,p,theta){dirac_trans <- function(n){dirac_bc_transprob(n,p)}
dirac_rate <- function(n){dirac_bc_rate(n,p)}
Lambda_seq_sim(Lambda_jc_sim(n,dirac_trans),theta,
               dirac_rate)}

#Simulate a Bolthausen-Sznitman n-coalescent with external branches prolongued by 1
#as appearing in the limit of Berestycki, Berestycki, Schweinsberg (Ann. Prob. 2013)
#
#Input conf1=Jump chain of an ordinary Bolthausen-Sznitman n-coalescent
#theta=scaled mutation rate, totratef=total rate function of 
#ordinary Bolthausen-Sznitman n-coalescent 
#Output: n times s matrix of zeros and ones, rows=individuals, cols=mutations 
bszplus1_seq_sim0 <- function(conf1,theta,totratef){
  n <- ncol(conf1) 
  seq1 <- NULL
  theta <- theta/2 #Mutations with half rate s.t. E(pw. diff. of 2 seqs)=\theta 
  bc_jc <- apply(conf1,1,function(x){length(unique(x))}) #Get block counts
  for (i in 1:(nrow(conf1)-1)){
    blocks <- unique(conf1[i,]) #Blocks present at jump
    n1 <- length(blocks) #number of blocks
    tot_rate <- totratef(n1) #total rate with n1 blocks
    wt <- rexp(1,tot_rate)
    mut1 <- rpois(n1,theta*wt)
    if (i == 1){mut1 <- mut1 + rpois(n1,theta)} #Numbers of mut's on the 
    #+1 branch parts are independent Poissonians
    for (j in 1:n1){
      if (mut1[j]>0){
        temp <- rep(0,n)
        temp[conf1[i,]==blocks[j]] <- 1
        temp <- matrix(temp,ncol=mut1[j],nrow=n)
        seq1 <- cbind(seq1,temp)}
    }
  }                 
  return(seq1)}
#Wrapper for producing sample of n-sequences under the model BSZ+1 given n, theta
bszplus1_seq_sim <- function(n,theta){
  bszplus1_seq_sim0(Lambda_jc_sim(n,bsz_bc_transprob),theta,
                    bsz_bc_rate)} 