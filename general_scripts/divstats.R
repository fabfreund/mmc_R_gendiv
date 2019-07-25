#' General input for all functions: SNP matrix 'seq1' of 0-1 (rows sequences, columns SNPs), 
#' quantile vector 'quant_v' 
#' Needs ape

#' Output: Quantiles 'quant_v' of allele frequencies
allele_freqs <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){if (!any(is.na(seq1))){
  if (!(is.matrix(seq1))){as.matrix(seq1)}  
  temp1 <- colSums(seq1)/nrow(seq1)
  temp2 <- quantile(temp1,quant_v)
  names(temp2) <- paste0("AF_qu",quant_v)
  return(temp2)
}
}


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
phylolength <- function(seq1,quant_v = c(.1,.3,.5,.7,.9)){
  if (is.matrix(seq1)){
    dist1 <- dist(seq1,method="manhattan")
    tree1 <- nj(dist1)  
    el1 <- tree1$edge.length} else {el1 <- 0}
  phy_l <- quantile(el1,quant_v)
  names(phy_l) <- paste0("phybl_qu",quant_v)
  return(phy_l)
}  

  
  
