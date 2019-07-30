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



#' Replicated from B. Eldon's work:
## Copyright (2014):  Bjarki Eldon 
## Distributed under GNU GPL v3 or later
## compile  .c file in terminal as: #c file in this folder, same copyright
## R CMD SHLIB -O3 CREBiepg.c -lm -lgsl -lgslcblas 
## code based on Polanski and Kimmel (2003)
## and Polanski, Bobrowski, Kimmel (2003)
## Nl is number leaves = sample size
## bp is the growth parameter (b)
## returns expected branch lengths E_b[B_i] where B_i
## is the random length of branches subtending i leaves

"REBiepg" <- function( Nl, bp )
{
  if( is.loaded( "CREBiepg.so" ) ){
    dyn.unload( "CREBiepg.so")
    dyn.load( "CREBiepg.so") }
  else
    dyn.load("CREBiepg.so")
  
  vx <- as.double( (1:(Nl)) * 0 )
  
  out =  .C("EBi_W", Nleaves = as.integer(Nl), betaparam = as.double(bp),   x = vx )
  
  return( out$x[-1] )
  
}
