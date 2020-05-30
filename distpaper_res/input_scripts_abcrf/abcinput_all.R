#' Use on sims with stats O, Ham, Phy, r^2, AF* (all with quants .1,.2,...,.9)
#' and full SFS
#' For abbreviations of stats, see Methods in our preprint

#' Third input in abcrf_xfold_2rep.R determines the stat sets used for ABC
#' Options are "af", "ohamaf", "fold", "fine", "sfs", "lda"
#' Set change_name to TRUE to add chosen option to filename  
#' see below for details

stats_used <- args1[3]
change_name <- TRUE


#' Rows in sims that are kept for each set of statistics used for model selection
o_coarse <- c(1,2,4,6,8,10:12)
o_fine <- 1:12
ham_coarse <- seq(13,21,2)
ham_fine <- 13:21
phy_coarse <- seq(22,30,2)
phy_fine <- 22:30
r2_coarse <- seq(31,39,2)  
r2_fine <- 31:39
s_pi <- c(40,41)
D_H <- c(42,43)
af_coarse <- seq(44,52,2)
af_fine <- 44:52
SFS <- c(53:nrow(sims1))

#' Different choices of sets of stats for the ABC, including transformations 
if (stats_used=="af"){
rows_abc <- list(c(s_pi,af_coarse),
                 c(o_coarse,s_pi,af_coarse),
                 c(ham_coarse,s_pi,af_coarse),
                 c(phy_coarse,s_pi,af_coarse),
                 c(r2_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,phy_coarse,r2_coarse),
                 c(ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse),
                 c(o_coarse,ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse))
#' Input vector
stats1 <- vector("list",9)
#' Which stat set is used for importance scoring?
index_imp <- 9
#' Names of sets of stats
names(stats1) <- c("AF+","AF+, O","AF+, Ham","AF+, Phy","AF+, r2",
                   "AF+,O,Ham","FULL - AF+","FULL - O","FULL")
}


if (stats_used=="ohamsfs"){
  sims2 <- apply(sims1[SFS[-(1:14)],],2,sum) 
  sims1 <- rbind(sims1,S15plus=sims2)
  SFSl <- c(SFS[1:14],max(SFS)+1)
  o_coarse <- o_coarse[2:6]
  rows_abc <- list(o_coarse,
                   ham_coarse,
                   SFSl,
                   c(o_coarse,ham_coarse),
                   c(o_coarse,SFSl),
                   c(ham_coarse,SFSl),
                   c(o_coarse,ham_coarse,SFSl))
  #' Input vector
  stats1 <- vector("list",7)
  #' Which stat set is used for importance scoring?
  index_imp <- 7
  #' Names of sets of stats
  names(stats1) <- c("O","Ham","SFSl",
                     "O, Ham","O,SFSl","Ham,SFSl","FULL")
}



#' Only use this if sample size n>30, since we lump 15+
#' 
if (stats_used=="fold"){
#' Function to fold SFS
  fold_f <- function(sfs1){l_sfs <- length(sfs1)
                           sfs1 <- sfs1 + sfs1[length(sfs1):1]
                           last_entry <- ceiling(l_sfs/2)
                           sfs1 <- sfs1[1:last_entry]
                           if (l_sfs %% 2 >0){
                    sfs1[last_entry] <- sfs1[last_entry]/2}
                           names(sfs1) <- paste0("fS",1:length(sfs1))
                           return(sfs1)}
#' Add folded SFS as further rows  
  sims2 <- apply(sims1[SFS,],2,fold_f)
  sims1 <- rbind(sims1,sims2)
#' Add lumped fSFS to sims1 
  sims3 <- apply(sims2[-(1:14),],2,sum) 
  sims1 <- rbind(sims1,fS15plus=sims3)
#' Write down rows in sims1 where the folded SFS is located
  fSFS <- max(SFS) + 1:nrow(sims2)
  fSFSl <- c(max(SFS) + c(1:14),max(fSFS)+1)
  Si_t <- c(max(SFS),max(fSFS)+1) 
  rows_abc <- list(c(s_pi,fSFS),
                   c(s_pi,fSFSl),
                   c(s_pi,Si_t),
                   c(ham_coarse,phy_coarse,r2_coarse),
                   c(ham_coarse,s_pi,Si_t),
                   c(phy_coarse,s_pi,Si_t),
                   c(r2_coarse,s_pi,Si_t),
                   c(ham_coarse,r2_coarse,s_pi,Si_t),
                   c(ham_coarse,phy_coarse,r2_coarse,s_pi,fSFS)
                   )
  #' Input vector
  stats1 <- vector("list",9)
  #' Which stat set is used for importance scoring?
  index_imp <- 9
  #' Names of sets of stats
  names(stats1) <- c("fSFS+","fSFSl+","fSi.t+",
                     "Ham,Phy,r2","Ham,fSi.t+",
                     "Phy,fSi.t+","r2,fSi.t+","H,r2,fSi.t+",
                     "FULL")
}


if (stats_used=="fine"){
  rows_abc <- list(c(ham_fine,phy_fine,r2_fine),
                   o_fine,c(s_pi,af_fine),
                   c(o_fine,ham_fine,s_pi,af_fine),
                   c(ham_fine,phy_fine,r2_fine,s_pi,af_fine),
                   c(o_fine,ham_fine,phy_fine,r2_fine),
                   c(o_fine,phy_fine,r2_fine,s_pi,af_fine),
                   c(o_fine,ham_fine,phy_fine,r2_fine,s_pi,af_fine))
  #' Input vector
  stats1 <- vector("list",8)
  #' Which stat set is used for importance scoring?
  index_imp <- 8
  #' Names of sets of stats
  names(stats1) <- c("Ham, Phy, r^2","O","AF+",
                     "O, Ham,AF+","FULL - O",
                     "FULL - AF+","FULL - Ham","FULL")
}

if (stats_used=="fine2"){
  rows_abc <- list(c(o_coarse,ham_coarse,s_pi,af_fine),
                   c(o_coarse,ham_coarse,s_pi,af_coarse))
  #' Input vector
  stats1 <- vector("list",2)
  #' Which stat set is used for importance scoring?
  index_imp <- 1
  #' Names of sets of stats
  names(stats1) <- c("O,Ham,fineAF+","O,Ham,AF+")
}


if (stats_used=="sfs"){
  sims2 <- apply(sims1[SFS[-(1:14)],],2,sum) 
  sims1 <- rbind(sims1,S15plus=sims2)
  SFSl <- c(SFS[1:14],max(SFS)+1)
  Si_t <- SFSl[c(1,15)]
  #index_imp <- 8
  rows_abc <- list(c(s_pi,D_H,SFS),
                   c(s_pi,D_H,SFSl),
                   SFS,
                   Si_t,
                   c(o_coarse,Si_t),
                   c(o_coarse,s_pi,D_H,Si_t),
                   c(o_coarse,ham_coarse,SFSl),
                   c(o_coarse,ham_coarse,s_pi,D_H,SFSl),
                   c(o_coarse,ham_coarse,phy_coarse,r2_coarse,s_pi,D_H,SFS)
                   )
  #' Input vector
  stats1 <- vector("list",9)
  #' Which stat set is used for importance scoring?
  index_imp <- 8
  #' Names of sets of stats
  names(stats1) <- c("SFS*","SFSl*","SFS","Si.t","O,Si.t","O,Si.t*","O,Ham,SFSl",
                     "O,Ham,SFSl*","FULL")
}

#' adds linear discriminant
if (stats_used=="lda"){
  #' To avoid problems when scaling constant variables
  scale_const <- function(v){if (sd(v,na.rm = TRUE)>0){out1 <- scale(v)} else {
                             out1 <- scale(v,scale = FALSE)}
                             }
  lda1 <- TRUE
  sims1 <- t(apply(sims1,1,scale_const)) #scale variables to interpret l.d. loadings
  rows_abc <- list(c(s_pi,af_coarse),
                   c(o_coarse,s_pi,af_coarse),
                   c(ham_coarse,s_pi,af_coarse),
                   c(phy_coarse,s_pi,af_coarse),
                   c(r2_coarse,s_pi,af_coarse),
                   c(o_coarse,ham_coarse,s_pi,af_coarse),
                   c(o_coarse,ham_coarse,phy_coarse,r2_coarse),
                   c(ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse),
                   c(o_coarse,ham_coarse,phy_coarse,r2_coarse,s_pi,af_coarse))
  #' Input vector
  stats1 <- vector("list",9)
  #' Which stat set is used for importance scoring?
  index_imp <- 9
  #' Names of sets of stats
  names(stats1) <- c("AF+","AF+, O","AF+, Ham","AF+, Phy","AF+, r2",
                     "AF+,O,Ham","FULL - AF+","FULL - O","FULL")
}
#' How many ABC runs per set of stats per replication
nreps <- 10
#'Output objects: oob errors and variable importances
res_conf <- vector("list",length(stats1))
names(res_conf) <- names(stats1)
res_imp <- vector("list",nreps)

