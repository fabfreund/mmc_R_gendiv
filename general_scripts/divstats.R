#' General input for all functions: SNP matrix 'seq1' of 0-1 (rows sequences, columns SNPs), 
#' quantile vector 'quant_v' 
#' Site frequency spectrum 'sfs', i.e. vector of integers (and zeros) 


#' Output: Quantiles 'quant_v' of allele frequencies
allele_freqs <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){if (!any(is.na(seq1))){
  if (!(is.matrix(seq1))){as.matrix(seq1)}  
  temp1 <- colSums(seq1)/nrow(seq1)
  temp2 <- quantile(temp1,quant_v)
  names(temp2) <- paste0("AF_qu",quant_v)
  return(temp2)
}
}

#' Output: Complete (polarized) site frequency spectrum
#' WARNING: THIS FUNCTION ONLY SHOWS THE RIGHT RESULT IF
#' THE MATRIX IS A PROPER MATRIX THAT ONLY SHOWS SNPs, SO
#' NO MONOMORPHIC POSITIONS (i.e. NO COLUMN WITH ONLY 1s OR ONLY 0s)
#' 2nd WARNING: NO NA ALLOWED, ELSE RETURNS ONLY ZEROS 
spectrum01 <- function(seq1){if (!any(is.na(seq1))){
  if (!(is.matrix(seq1))){seq1 <- as.matrix(seq1)}  
  n1 <- nrow(seq1) #rows are individuals
  b <-rep(0,n1-1) #only segregating sites are simulated
  a <- colSums(seq1)
  #for each mutation, count it for the right multiplicity
  for (i in seq(along=a)) {b[a[i]]<-b[a[i]]+1}} else {b <- rep(0,n1-1)} 
  names(b) <- paste0("S",1:(n1-1)) 
  return(b)}

#' Further input! fold1: Fold SFS?; scale1: Scale sfs/sum(sfs)?; lump: Should the high classes
#' of the SFS be added together, i.e. SFS (S1,S2,...Sn-1) -> (S1,...,Sln-1,sum(Sln,...,Sn-1))
#' If fold1=TRUE and lump=TRUE, lumping is after folding
fnl_sfs <- function(sfs1,fold1=TRUE,scale1=TRUE,lump=FALSE,ln=2){
  s1 <- sum(sfs1);l <- length(sfs1)
  if (lump){
    if ((ln > l & !fold1) | (ln >  ceiling(l/2) & fold1)){warning("bad lumping");return(NA)}}
  if (scale1 & sum(s1 > 0)){sfs1 <- sfs1/s1}
  if (fold1){b <- rep(0,ceiling(l/2))
  if (s1>0){for (i in seq(along=b[-1])) {b[i] <- sfs1[i]+sfs1[l+1-i]}
    if (l/2 > floor(l/2)){b[ceiling(l/2)] <- sfs1[ceiling(l/2)]} else {b[l/2] <- sfs1[l/2]+sfs1[l/2+1]}}
  sfs1 <- b}
  if (lump){sfs1 <- c(sfs1[1:(ln-1)],sum(sfs1[-(1:(ln-1))]))}
  names(sfs1) <- paste0("S",1:length(sfs1))
  if (lump){
    names(sfs1) <- c(paste0("S",1:(length(sfs1)-1)),
                     paste0("S",length(sfs1),"+"))}
  return(sfs1)}

#' Output: Quantiles 'quant_v' of Hamming distances of all pairs of sequences
hammfun <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  if (is.matrix(seq1)){
    dist1 <- dist(seq1,method="manhattan")
  } else {dist1 <- 0}
  gendist <- quantile(dist1,quant_v)
  names(gendist) <- paste0("hammd_qu",quant_v)
  return(gendist)
}  

#' Output: Quantiles 'quant_v' of (Hamming distances)/(nucleotide diversity) of all pairs of sequences
hammfun_sc <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  if (is.matrix(seq1)){
    dist1 <- dist(seq1,method="manhattan")
    dist1 <- dist1/mean(dist1)   
  } else {dist1 <- 0}
  gendist <- quantile(dist1,quant_v)
  names(gendist) <- paste0("hammd_sc_qu",quant_v)
  return(gendist)
}  

#' Output: Quantiles 'quant_v' of branch lengths of a neighbour-joining reconstruction of the genealogy/phylogeny 
#' of 'seq1'
#' needs R package 'ape'
#' require(ape) not in function for speed issues
require(ape)
phylolength <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  if (is.matrix(seq1)){
    dist1 <- dist(seq1,method="manhattan")
    tree1 <- nj(dist1)  
    el1 <- tree1$edge.length} else {el1 <- 0}
  phy_l <- quantile(el1,quant_v)
  names(phy_l) <- paste0("phybl_qu",quant_v)
  return(phy_l)
}  


#' Output: vector (nucleotide diversity, Number segregating sites) 
f_nucdiv_S <- function(sfs){n1 <- length(sfs)+1
PI <- sum((1:(n1-1))*((n1-1):1)*sfs)
pi <- 2*PI/(n1*(n1-1))
names(pi) <- "nucdiv" 
S <- sum(sfs)
names(S) <- "S"
return(c(pi,S))}  

#' Output of r^2 as LD measure between pairs of SNPs
r2fun <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  cor1 <- cor(seq1)
  cor1 <- cor1[lower.tri(seq1)]
  r2 <- cor1^2
  return(c(r2=quantile(r2,quant_v,na.rm = TRUE)))
}

#' Internal function that computes the minimal observable clade size of the sequence
#' represented by row i
OCi <- function(i,seq1){
  targetseq <- seq1[i,]
  mut_pos <- which(targetseq==1) #All mutations of seq. i
  if (length(mut_pos)<1){return(nrow(seq1))} else {
    seq2 <- seq1[,mut_pos] #seq2 shows only mutation pos. of seq. i
    if (length(mut_pos)==1){return(sum(seq2))} else {
      return(min(colSums(seq2)))} #How many seqs have all these muts?
  }}


#' Output: Vector of all minimal observable clade sizes (of row 1,...,n)
vect_OCi <- function(seq1){if (!all(is.na(seq1))){
  if (!(is.matrix(seq1))){as.matrix(seq1)}  
  log1 <- (colSums(seq1)>1)
  if (sum(log1)==0){out1 <- rep(nrow(seq1),nrow(seq1))} else {
    seq1 <- seq1[,log1] #Remove all private mutations
    if (sum(log1)==1){seq1 <- as.matrix(seq1)}
    out1 <- sapply(1:nrow(seq1),OCi,seq1)}} else {out1 <- rep(n,n)}
  return(out1)}

#' Output: Harmonic mean and quantiles 'quant_v' of minimal observable clade sizes
#' uses 'psych' package for harmonic mean
require(psych)
quant_hm_oc <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  oc_vals <- vect_OCi(seq1)
  out1 <- c(harmonic.mean(oc_vals),quantile(oc_vals,quant_v))
  names(out1) <- c("oc_hm",paste0("oc_",quant_v))
  return(out1)
} 

#' Output: Sample mean and standard deviation of minimal observable clade sizes
mean_sd_oc <- function(seq1){
  oc_vals <- vect_OCi(seq1)
  out1 <- c(mean(oc_vals),sd(oc_vals))
  names(out1) <- c("oc_mean","oc_sd")
  return(out1)
} 

#' Output: Quantiles 'quant_v' of scaled observable minimal clade sizes O(i)/max_i(O(i))
quant_oc_sc <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  oc_vals <- vect_OCi(seq1)
  oc_vals <- oc_vals/max(oc_vals)
  out1 <- c(quantile(oc_vals,quant_v))
  names(out1) <- c(paste0("oc_sc_",quant_v))
  return(out1)
} 


