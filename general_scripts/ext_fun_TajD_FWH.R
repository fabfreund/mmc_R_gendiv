#' To compute Tajima's D and Fay and Wu's H, we load code from  
#' sdwfrost/popseq: Population Analysis of Deep Sequencing Data
#' available on GitHub https://github.com/sdwfrost/popseq

require(devtools)
source_url('https://raw.githubusercontent.com/sdwfrost/popseq/master/R/sfsR.R')

#' sfsR uses 'hapmatrix' as argument, which is like our 0-1 SNP matrix, but with
#' mutations coded as 2 and wild type as 1. We provide 
#' Needs R package e1071


require(e1071)

#' Input: 0-1 SNP matrix Output: named vector (TajimaD, FayandWuH)
D_H <- function(seq1){unlist(sfsR(seq1+1))}

