source("../general_scripts/divstats.R")
source("../general_scripts/ext_fun_TajD_FWH.R")

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


#' This wrapper needs the sample size n_ind as additional variable
divfun_most <- function(seq1,n_ind){
  quantv1 <- seq(.1,.9,.1)
  if (is.matrix(seq1)){
    out1 <- c(quant_hm_oc(seq1,quant_v = quantv1),mean_sd_oc(seq1),
              hammfun(seq1,quant_v = quantv1),phylolength(seq1,quant_v = quantv1),r2fun(seq1,quant_v = quantv1),
              f_nucdiv_S(spectrum01(seq1)),D_H(seq1),allele_freqs(seq1,quant_v = quantv1),spectrum01(seq1))}
  else {out1 <- rep(NA,51+n_ind)}
  names(out1) <-c("hm(O)",paste("O: qu",seq(.1,.9,.1)),"mean(O)","sd(O)",
                  paste("Ham: qu",seq(.1,.9,.1)),
                  paste("Phy: qu",seq(.1,.9,.1)),
                  paste("r2: qu",seq(.1,.9,.1)),"Nucl. div.","S","Tajima's D", "Fay & Wu's H",
                  paste("AF: qu",seq(.1,.9,.1)),paste0("S",1:(n_ind-1)))
  return(out1)}
  
  
  divfun_no_s <- function(seq1){
    if (is.matrix(seq1)){
      log1 <- (colSums(seq1)>1)
      if (sum(log1)==0){out1 <- rep(NA,25)} else {
        seq1 <- as.matrix(seq1[,log1])
        out1 <- c(quant_hm_oc(seq1),mean_sd_oc(seq1),
                   hammfun(seq1),r2fun(seq1),
                   f_nucdiv_S(spectrum01(seq1)),allele_freqs(seq1))}}
    else {out1 <- rep(NA,25)}
    names(out1) <-c("hm(O)",paste("O: qu",seq(.1,.9,.2)),"mean(O)","sd(O)",
                    paste("Ham: qu",seq(.1,.9,.2)),
                    paste("r2: qu",seq(.1,.9,.2)),
                    "Nucl. div.","S",
                    paste("AF: qu",seq(.1,.9,.2)))
    return(out1)}
  
  divfun_robust <- function(seq1,n_ind){
    if (is.matrix(seq1)){
      out1 <- c(quant_oc_sc(seq1),
                hammfun_sc(seq1),
                fnl_sfs(spectrum_fab(seq1),FALSE,TRUE,TRUE,15))
    } else {out1 <- rep(NA,25)}
    names(out1) <-c("oc_sc_.1","oc_sc_.3","oc_sc_.5","oc_sc_.7","oc_sc_.9",
                    "ham_sc_.1","ham_sc_.3","ham_sc_.5","ham_sc_.7","ham_sc_.9", 
                    paste0("scS",1:14),"scS15+")
    
    return(out1)}
  
  
  